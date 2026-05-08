struct SmithConfig
    max_probability::Float64
    min_probability::Float64
end

SmithConfig() = SmithConfig(1.0, 0.0)

mutable struct SmithPacer
    config::SmithConfig
    last_time::Float64
    last_probability::Float64
    delayed_budget::Float64
    delayed_spent::Float64
    delay_buffer::Vector{Tuple{Float64,Float64,Float64}}  # (time, budget, spent)
    τ::Float64
end

SmithPacer(now::Float64; τ::Float64=5.0) =
    SmithPacer(SmithConfig(), now, 0.0, 0.0, 0.0, Tuple{Float64,Float64,Float64}[], τ)

struct SmithResult
    probability::Float64
    predicted_balance::Float64
    ideal_rate::Float64
    time_remaining::Float64
end

function calculate_smith_mode!(pacer::SmithPacer;
        current_time::Float64,
        start_time::Float64,
        end_time::Float64,
        total_budget::Float64,
        current_spent::Float64,
        max_spend_rate::Float64)::SmithResult

    dt = current_time - pacer.last_time

    # Guard: not in flight
    if current_time < start_time || current_time >= end_time
        return SmithResult(0.0, 0.0, 0.0, 0.0)
    end

    # Guard: no eligible spend
    if max_spend_rate <= 0.0
        return SmithResult(0.0, 0.0, 0.0, 0.0)
    end

    # Guard: dt too small
    if dt < 0.5
        return SmithResult(pacer.last_probability, 0.0, 0.0, 0.0)
    end

    current_balance = total_budget - current_spent

    # Record current observation in delay buffer
    push!(pacer.delay_buffer, (current_time, current_balance, current_spent))

    # Find delayed observation (τ seconds ago)
    delay_target = current_time - pacer.τ
    delayed_balance = current_balance
    delayed_spent = current_spent

    for i in length(pacer.delay_buffer):-1:1
        t_i, b_i, s_i = pacer.delay_buffer[i]
        if t_i <= delay_target
            delayed_balance = b_i
            delayed_spent = s_i
            break
        end
    end

    # Prune old entries beyond 2τ
    cutoff = current_time - 2.0 * pacer.τ
    filter!(entry -> entry[1] >= cutoff, pacer.delay_buffer)

    # Smith predictor: reconstruct current balance from delayed observation
    # B̂(t) = B(t-τ) - (S(t) - S(t-τ))
    predicted_balance = max(0.0, delayed_balance - (current_spent - delayed_spent))

    time_remaining = end_time - current_time
    if time_remaining <= 0.0
        return SmithResult(0.0, predicted_balance, 0.0, 0.0)
    end

    # Ideal spend rate to exhaust remaining balance over remaining time
    ideal_rate = predicted_balance / time_remaining

    # Probability = idealRate / maxSpendRate
    # "What fraction of available supply should we accept?"
    probability = clamp(ideal_rate / max_spend_rate, 0.0, 1.0)

    # Update state
    pacer.last_time = current_time
    pacer.last_probability = probability

    return SmithResult(probability, predicted_balance, ideal_rate, time_remaining)
end
