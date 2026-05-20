struct MPCConfig
    max_probability::Float64
    min_probability::Float64
end

MPCConfig() = MPCConfig(1.0, 0.0)

mutable struct MPCPacer
    config::MPCConfig
    last_time::Float64
    last_probability::Float64
    delayed_budget::Float64
    delayed_spent::Float64
    delay_buffer::Vector{Tuple{Float64,Float64,Float64}}  # (time, budget, spent)
    τ::Float64
    α_hat::Float64
end

MPCPacer(now::Float64; τ::Float64=5.0) =
    MPCPacer(MPCConfig(), now, 0.0, 0.0, 0.0, Tuple{Float64,Float64,Float64}[], τ, 1.0)

struct MPCResult
    probability::Float64
    predicted_balance::Float64
    observed_balance::Float64
    ideal_rate::Float64
    safe_rate::Float64
    α_hat::Float64
    time_remaining::Float64
end

function calculate_mpc_mode!(pacer::MPCPacer;
        current_time::Float64,
        start_time::Float64,
        end_time::Float64,
        total_budget::Float64,
        current_spent::Float64,
        max_spend_rate::Float64)::MPCResult

    dt = current_time - pacer.last_time

    # Guard: not in flight
    if current_time < start_time || current_time >= end_time
        return MPCResult(0.0, 0.0, 0.0, 0.0, 0.0, pacer.α_hat, 0.0)
    end

    # Guard: no eligible spend
    if max_spend_rate <= 0.0
        return MPCResult(0.0, 0.0, 0.0, 0.0, 0.0, pacer.α_hat, 0.0)
    end

    # Guard: dt too small
    if dt < 0.5
        return MPCResult(pacer.last_probability, 0.0, 0.0, 0.0, 0.0, pacer.α_hat, 0.0)
    end

    current_balance = total_budget - current_spent

    # Record current observation in delay buffer
    push!(pacer.delay_buffer, (current_time, current_balance, current_spent))

    # Find delayed observation (τ seconds ago) and double-delayed (2τ ago)
    delay_target = current_time - pacer.τ
    delay_target_2 = current_time - 2.0 * pacer.τ

    delayed_balance = current_balance
    delayed_spent = current_spent
    double_delayed_balance = current_balance
    double_delayed_spent = current_spent

    for i in length(pacer.delay_buffer):-1:1
        t_i, b_i, s_i = pacer.delay_buffer[i]
        if t_i <= delay_target
            delayed_balance = b_i
            delayed_spent = s_i
            break
        end
    end

    for i in length(pacer.delay_buffer):-1:1
        t_i, b_i, s_i = pacer.delay_buffer[i]
        if t_i <= delay_target_2
            double_delayed_balance = b_i
            double_delayed_spent = s_i
            break
        end
    end

    # Prune old entries beyond 3τ
    cutoff = current_time - 3.0 * pacer.τ
    filter!(entry -> entry[1] >= cutoff, pacer.delay_buffer)

    # Estimate α̂: ratio of observed drain to predicted drain
    # observed drain = how much balance actually dropped between 2τ and τ ago
    # predicted drain = how much we granted (spent) in that same interval
    predicted_drain = delayed_spent - double_delayed_spent
    observed_drain = double_delayed_balance - delayed_balance
    if predicted_drain > 1e-6
        pacer.α_hat = clamp(observed_drain / predicted_drain, 0.01, 2.0)
    end

    # Adaptive Smith estimate: predict current balance
    # B̂(t) = B(t-τ) - α̂ · (S(t) - S(t-τ))
    predicted_balance = max(0.0, delayed_balance - pacer.α_hat * (current_spent - delayed_spent))

    # Observed balance from delayed signal (safe/naive estimate)
    observed_balance = max(0.0, delayed_balance)

    time_remaining = end_time - current_time
    if time_remaining <= 0.0
        return MPCResult(0.0, predicted_balance, observed_balance, 0.0, 0.0, pacer.α_hat, 0.0)
    end

    # MPC: unconstrained rate from adaptive Smith, capped by safe naive rate
    rate_free = predicted_balance / time_remaining
    rate_safe = observed_balance / time_remaining
    ideal_rate = min(rate_free, rate_safe)

    # Probability = idealRate / maxSpendRate
    probability = clamp(ideal_rate / max_spend_rate, 0.0, 1.0)

    # Update state
    pacer.last_time = current_time
    pacer.last_probability = probability

    return MPCResult(probability, predicted_balance, observed_balance, rate_free, rate_safe, pacer.α_hat, time_remaining)
end
