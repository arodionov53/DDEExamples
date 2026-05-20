struct ArtsteinConfig
    max_probability::Float64
    min_probability::Float64
end

ArtsteinConfig() = ArtsteinConfig(1.0, 0.0)

mutable struct ArtsteinPacer
    config::ArtsteinConfig
    last_time::Float64
    last_probability::Float64
    τ::Float64
    # Full history buffer: (time, balance, spent, rate)
    # rate = instantaneous spend rate at that tick
    history::Vector{Tuple{Float64,Float64,Float64,Float64}}
    cumulative_integral::Float64
end

ArtsteinPacer(now::Float64; τ::Float64=5.0) =
    ArtsteinPacer(ArtsteinConfig(), now, 0.0, τ,
                  Tuple{Float64,Float64,Float64,Float64}[], 0.0)

struct ArtsteinResult
    probability::Float64
    predicted_balance::Float64
    integral_correction::Float64
    ideal_rate::Float64
    time_remaining::Float64
end

function calculate_artstein_mode!(pacer::ArtsteinPacer;
        current_time::Float64,
        start_time::Float64,
        end_time::Float64,
        total_budget::Float64,
        current_spent::Float64,
        max_spend_rate::Float64)::ArtsteinResult

    dt = current_time - pacer.last_time

    # Guard: not in flight
    if current_time < start_time || current_time >= end_time
        return ArtsteinResult(0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Guard: no eligible spend
    if max_spend_rate <= 0.0
        return ArtsteinResult(0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Guard: dt too small
    if dt < 0.5
        return ArtsteinResult(pacer.last_probability, 0.0, 0.0, 0.0, 0.0)
    end

    current_balance = total_budget - current_spent

    # Compute instantaneous spend rate from last observation
    instantaneous_rate = 0.0
    if !isempty(pacer.history)
        prev_t, _, prev_spent, _ = pacer.history[end]
        dt_obs = current_time - prev_t
        if dt_obs > 0.0
            instantaneous_rate = (current_spent - prev_spent) / dt_obs
        end
    end

    # Record current observation
    push!(pacer.history, (current_time, current_balance, current_spent, instantaneous_rate))

    # Find delayed balance observation B(t-τ)
    delay_target = current_time - pacer.τ
    delayed_balance = current_balance

    for i in length(pacer.history):-1:1
        t_i, b_i, _, _ = pacer.history[i]
        if t_i <= delay_target
            delayed_balance = b_i
            break
        end
    end

    # Artstein integral: ∫_{t-τ}^{t} rate(s) ds
    # Computed via trapezoidal rule over the history buffer within [t-τ, t]
    integral_correction = 0.0
    for i in 2:length(pacer.history)
        t_prev, _, _, rate_prev = pacer.history[i-1]
        t_curr, _, _, rate_curr = pacer.history[i]

        # Clip interval to [delay_target, current_time]
        t0 = max(t_prev, delay_target)
        t1 = min(t_curr, current_time)
        if t1 > t0
            # Trapezoidal rule
            integral_correction += 0.5 * (rate_prev + rate_curr) * (t1 - t0)
        end
    end

    # Artstein prediction: B̂(t) = B(t-τ) - ∫_{t-τ}^{t} rate(s) ds
    predicted_balance = max(0.0, delayed_balance - integral_correction)

    # Prune old entries beyond 2τ
    cutoff = current_time - 2.0 * pacer.τ
    filter!(entry -> entry[1] >= cutoff, pacer.history)

    time_remaining = end_time - current_time
    if time_remaining <= 0.0
        return ArtsteinResult(0.0, predicted_balance, integral_correction, 0.0, 0.0)
    end

    # Ideal spend rate to exhaust predicted remaining balance
    ideal_rate = predicted_balance / time_remaining

    # Probability = idealRate / maxSpendRate
    probability = clamp(ideal_rate / max_spend_rate, 0.0, 1.0)

    # Update state
    pacer.last_time = current_time
    pacer.last_probability = probability

    return ArtsteinResult(probability, predicted_balance, integral_correction, ideal_rate, time_remaining)
end
