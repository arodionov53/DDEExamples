struct CorrDenomConfig
    max_probability::Float64
    min_probability::Float64
end

CorrDenomConfig() = CorrDenomConfig(1.0, 0.0)

mutable struct CorrDenomPacer
    config::CorrDenomConfig
    last_time::Float64
    last_probability::Float64
    τ::Float64
end

CorrDenomPacer(now::Float64; τ::Float64=5.0) =
    CorrDenomPacer(CorrDenomConfig(), now, 0.0, τ)

struct CorrDenomResult
    probability::Float64
    observed_balance::Float64
    ideal_rate::Float64
    time_remaining::Float64
end

function calculate_corrdenom_mode!(pacer::CorrDenomPacer;
        current_time::Float64,
        start_time::Float64,
        end_time::Float64,
        total_budget::Float64,
        current_spent::Float64,
        max_spend_rate::Float64)::CorrDenomResult

    dt = current_time - pacer.last_time

    # Guard: not in flight
    if current_time < start_time || current_time >= end_time
        return CorrDenomResult(0.0, 0.0, 0.0, 0.0)
    end

    # Guard: no eligible spend
    if max_spend_rate <= 0.0
        return CorrDenomResult(0.0, 0.0, 0.0, 0.0)
    end

    # Guard: dt too small
    if dt < 0.5
        return CorrDenomResult(pacer.last_probability, 0.0, 0.0, 0.0)
    end

    observed_balance = total_budget - current_spent

    # Corrected denominator: divide observed balance by (time_remaining + τ)
    # This compensates for the fact that the observation is τ seconds stale
    time_remaining = end_time - current_time
    if time_remaining <= 0.0
        return CorrDenomResult(0.0, observed_balance, 0.0, 0.0)
    end

    ideal_rate = observed_balance / (time_remaining + pacer.τ)

    # Probability = idealRate / maxSpendRate
    probability = clamp(ideal_rate / max_spend_rate, 0.0, 1.0)

    # Update state
    pacer.last_time = current_time
    pacer.last_probability = probability

    return CorrDenomResult(probability, observed_balance, ideal_rate, time_remaining)
end
