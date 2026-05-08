struct CampaignMetadata
    total_budget::Int64
    start_time::Float64
    end_time::Float64
    avg_cost_per_impression::Int64
end

mutable struct CampaignState
    current_time::Float64
    latest_eligible_impression_rate::Int64
    last_observed_eligible_impression_rate::Int64
    win_percent::Float64
    throttle::Float64
    spent_impressions::Int64
    spent_budget::Int64
    pending_wins::PendingWinQueue
    pending_rates::PendingRateQueue
end

struct SimulationExpects
    ignore_pacing_violations_while_in_flight::Bool
    pacing_delta_percentage::Float64
    expected_impressions::Union{Nothing,Int64}
    expected_impressions_delta_percentage::Float64
end

struct SimulationConfig
    estimated_eligible_impression_rate::Int64
    estimated_eligible_impression_rate_stddev::Int64
    estimated_eligible_impression_spike_percent::Float64
    eligible_impressions_generator::Function
    estimated_win_percent::Float64
    estimated_win_percent_stddev::Float64
    win_percent_generator::Function
    estimated_tick_interval::Float64
    estimated_tick_interval_stddev::Float64
    tick_interval_generator::Function
    tick_count::Int64
    tick_until_campaign_ends::Bool
    win_delay_generator::Union{Nothing,Function}
    eligible_impression_latency_generator::Union{Nothing,Function}
end

mutable struct SimulationUseCase
    name::String
    metadata::CampaignMetadata
    state::CampaignState
    sim_config::SimulationConfig
    expects::SimulationExpects
end

struct PlotDataPoint
    time::Float64
    tick_delta::Float64
    linear_on_target_ratio::Float64
    expected_budget_utilization::Float64
    budget_utilization::Float64
    throttle::Float64
    eligible_impression_rate::Int64
end

# ── SimulationUseCase Helper Methods ─────────────────────────────────────────

function elapsed_duration_percent(uc::SimulationUseCase)::Float64
    if uc.state.current_time < uc.metadata.start_time
        return 0.0
    end
    if uc.state.current_time >= uc.metadata.end_time
        return 1.0
    end
    total = uc.metadata.end_time - uc.metadata.start_time
    total <= 0.0 && return NaN
    elapsed = uc.state.current_time - uc.metadata.start_time
    return elapsed / total
end

function target_spent_budget(uc::SimulationUseCase)::Int64
    round(Int64, elapsed_duration_percent(uc) * Float64(uc.metadata.total_budget))
end

function eligible_spend_budget_rate(uc::SimulationUseCase)::Int64
    uc.state.last_observed_eligible_impression_rate * uc.metadata.avg_cost_per_impression
end

function projected_outstanding_spent_impression_rate(uc::SimulationUseCase)::Int64
    round(Int64, uc.state.throttle * Float64(uc.state.last_observed_eligible_impression_rate))
end

function projected_outstanding_spent_budget_rate(uc::SimulationUseCase)::Int64
    projected_outstanding_spent_impression_rate(uc) * uc.metadata.avg_cost_per_impression
end

function total_spent_exposure_rate(uc::SimulationUseCase)::Int64
    pending = total_budget(uc.state.pending_wins)
    uc.state.spent_budget + pending + projected_outstanding_spent_budget_rate(uc)
end

function actualized_spent_budget_percent(uc::SimulationUseCase)::Float64
    Float64(uc.state.spent_budget) / Float64(uc.metadata.total_budget)
end

function expected_budget_utilization(uc::SimulationUseCase)::Float64
    elapsed_duration_percent(uc)
end

function linear_on_target_budget_utilization(uc::SimulationUseCase)::Float64
    edp = elapsed_duration_percent(uc)
    edp <= 0.0 && return NaN
    actualized_spent_budget_percent(uc) / edp
end

function should_tick(uc::SimulationUseCase, i::Int)::Bool
    if uc.sim_config.tick_until_campaign_ends
        if uc.state.current_time < uc.metadata.end_time
            return true
        end
        campaign_ended = uc.state.current_time >= uc.metadata.end_time
        queue_empty = length(uc.state.pending_wins) == 0
        if queue_empty && campaign_ended
            return false
        end
        return true
    end
    return i < uc.sim_config.tick_count
end

# ── Simulation Tick Engine ───────────────────────────────────────────────────

function tick_with!(uc::SimulationUseCase, interval::Float64)
    # 1. Advance clock
    uc.state.current_time += interval

    now = uc.state.current_time
    has_delayed_feedback = !isnothing(uc.sim_config.win_delay_generator)

    # 2. Deliver pending wins that have arrived
    if has_delayed_feedback
        arrived = dequeue_until!(uc.state.pending_wins, now)
        for win in arrived
            uc.state.spent_impressions += win.impressions
            uc.state.spent_budget += win.budget
        end
    end

    # 3. Generate stochastic market conditions
    has_rate_latency = !isnothing(uc.sim_config.eligible_impression_latency_generator)

    actual_rate = uc.sim_config.eligible_impressions_generator(now, uc.sim_config)
    uc.state.latest_eligible_impression_rate = actual_rate

    if has_rate_latency
        latency = uc.sim_config.eligible_impression_latency_generator(now, uc.sim_config)
        obs = RateObservation(actual_rate, now, now + latency)
        enqueue!(uc.state.pending_rates, obs)
        rate, ok = deliver_latest!(uc.state.pending_rates, now)
        if ok
            uc.state.last_observed_eligible_impression_rate = rate
        end
    else
        uc.state.last_observed_eligible_impression_rate = actual_rate
    end

    uc.state.win_percent = uc.sim_config.win_percent_generator(now, uc.sim_config)

    # 4. Create new win events (using ground truth rate)
    actual_projected_imps = round(Int64, uc.state.throttle * Float64(actual_rate))
    actual_projected_budget = actual_projected_imps * uc.metadata.avg_cost_per_impression
    new_impressions = round(Int64, Float64(actual_projected_imps) * interval * uc.state.win_percent)
    new_budget = round(Int64, Float64(actual_projected_budget) * interval * uc.state.win_percent)

    if has_delayed_feedback && new_impressions > 0
        # Scale batches with interval
        interval_minutes = interval / 60.0
        num_batches = clamp(round(Int64, 3 + sqrt(interval_minutes * 4)), 3, 20)
        num_batches = min(num_batches, new_impressions)
        step_imps = div(new_impressions, num_batches)
        step_budget = div(new_budget, num_batches)

        budget_used = Int64(0)
        i = Int64(0)
        while i < new_impressions
            batch_imps = step_imps
            if i + step_imps > new_impressions
                batch_imps = new_impressions - i
            end

            batch_budget = step_budget
            if i + step_imps > new_impressions
                batch_budget = new_budget - budget_used
            end

            delay = uc.sim_config.win_delay_generator(now, uc.sim_config)
            event = SpendEvent(batch_imps, batch_budget, now, now + delay)
            enqueue!(uc.state.pending_wins, event)
            budget_used += batch_budget
            i += step_imps
        end
    else
        uc.state.spent_impressions += new_impressions
        uc.state.spent_budget += new_budget
    end
end

function tick!(uc::SimulationUseCase)::Float64
    interval = uc.sim_config.tick_interval_generator(uc.state.current_time, uc.sim_config)
    tick_with!(uc, interval)
    return interval
end

function update_throttle!(uc::SimulationUseCase, throttle::Float64)
    uc.state.throttle = throttle
end

include("check_constraints.jl")
