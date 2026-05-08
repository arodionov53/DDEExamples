struct PIDConfig
    kp::Float64
    ki::Float64
    kd::Float64
    max_probability::Float64
    min_probability::Float64
end

PIDConfig(; kp=1.0, ki=0.3, kd=0.05) = PIDConfig(kp, ki, kd, 1.0, 0.0)

mutable struct PIDPacer
    config::PIDConfig
    integral::Float64
    last_error::Float64
    last_time::Float64
    last_probability::Float64
end

PIDPacer(now::Float64; kp=1.0, ki=0.3, kd=0.05) =
    PIDPacer(PIDConfig(; kp, ki, kd), 0.0, 0.0, now, 0.0)

struct PIDResult
    probability::Float64
    error::Float64
    proportional_term::Float64
    integral_term::Float64
    derivative_term::Float64
    control_variable::Float64
end

function calculate_cruise_mode!(pacer::PIDPacer;
        current_time::Float64,
        start_time::Float64,
        end_time::Float64,
        target_spent::Float64,
        total_exposure::Float64,
        max_spend_rate::Float64)::PIDResult

    dt = current_time - pacer.last_time

    # Guard: not in flight
    if current_time < start_time || current_time >= end_time
        return PIDResult(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    # Guard: no eligible spend
    if max_spend_rate <= 0.0
        return PIDResult(0.0, pacer.last_error, 0.0, pacer.config.ki * pacer.integral, 0.0, 0.0)
    end

    # Guard: dt too small (< 0.5 seconds)
    if dt < 0.5
        return PIDResult(pacer.last_probability, pacer.last_error, 0.0, pacer.config.ki * pacer.integral, 0.0, 0.0)
    end

    # Error: level-based (cumulative dollars)
    error_val = target_spent - total_exposure

    # Proportional term
    p_term = pacer.config.kp * error_val

    # Integral with windup protection
    pacer.integral += error_val * dt
    max_integral = 10.0 / max(pacer.config.ki, 0.01)
    pacer.integral = clamp(pacer.integral, -max_integral, max_integral)
    i_term = pacer.config.ki * pacer.integral

    # Derivative term
    derivative = (error_val - pacer.last_error) / dt
    d_term = pacer.config.kd * derivative

    # Combined control variable
    cv = p_term + i_term + d_term

    # Output normalization → probability
    probability = clamp(cv / (max_spend_rate * dt), 0.0, 1.0)

    # Update state
    pacer.last_error = error_val
    pacer.last_time = current_time
    pacer.last_probability = probability

    return PIDResult(probability, error_val, p_term, i_term, d_term, cv)
end
