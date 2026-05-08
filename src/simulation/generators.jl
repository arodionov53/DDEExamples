function default_eligible_impression_generator(now::Float64, config::SimulationConfig)::Int64
    mean_val = Float64(config.estimated_eligible_impression_rate)
    stddev = Float64(config.estimated_eligible_impression_rate_stddev)
    result = max(Int64(0), round(Int64, randn() * stddev + mean_val))
    if rand() < config.estimated_eligible_impression_spike_percent
        return Int64(200) * result
    end
    return result
end

function spike_and_trickle_generator(now::Float64, config::SimulationConfig)::Int64
    mean_val = Float64(config.estimated_eligible_impression_rate)
    stddev = Float64(config.estimated_eligible_impression_rate_stddev)
    if rand() < config.estimated_eligible_impression_spike_percent
        return Int64(200) * max(Int64(0), round(Int64, randn() * stddev + mean_val))
    end
    return round(Int64, mean_val)
end

function default_win_rate_generator(now::Float64, config::SimulationConfig)::Float64
    clamp(randn() * config.estimated_win_percent_stddev + config.estimated_win_percent, 0.0, 1.0)
end

function default_tick_interval_generator(now::Float64, config::SimulationConfig)::Float64
    max(0.05, randn() * config.estimated_tick_interval_stddev + config.estimated_tick_interval)
end

function default_eligible_impression_latency_generator(now::Float64, config::SimulationConfig)::Float64
    clamp(exp(randn() * 0.6 + 1.6), 30.0, 300.0)
end

function load_event_latency_generator(now::Float64, config::SimulationConfig)::Float64
    -log(rand()) / 1.85
end

function impression_event_latency_generator(now::Float64, config::SimulationConfig)::Float64
    if rand() < 0.48
        exp(randn() * 1.1 + 2.9)
    else
        exp(randn() * 1.2 + 5.6)
    end
end

function tracker_event_latency_generator(now::Float64, config::SimulationConfig)::Float64
    exp(randn() * 1.08 + 5.45)
end

function click_event_latency_generator(now::Float64, config::SimulationConfig)::Float64
    r = rand()
    if r < 0.45
        exp(randn() * 1.5 + 5.2)
    elseif r < 0.70
        exp(randn() * 0.9 + 8.4)
    else
        exp(randn() * 0.8 + 10.1)
    end
end
