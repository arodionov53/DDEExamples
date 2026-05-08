function check_hard_constraints(uc::SimulationUseCase, result::PIDResult)::Tuple{Bool,String}
    # Check 1: Probability zero outside flight
    if uc.state.current_time < uc.metadata.start_time ||
       uc.state.current_time >= uc.metadata.end_time
        if uc.state.throttle != 0.0
            return (false, "throttle $(uc.state.throttle) should be 0 outside campaign flight")
        end
    end

    # Check 2: Probability bounds
    if !isnan(uc.state.throttle) && (uc.state.throttle < 0.0 || uc.state.throttle > 1.0)
        return (false, "throttle $(uc.state.throttle) out of range [0.0, 1.0]")
    end

    # Check 3: Budget constraint (5% tolerance)
    budget_with_tolerance = Float64(uc.metadata.total_budget) * 1.05
    if Float64(uc.state.spent_budget) > budget_with_tolerance
        overspend = Float64(uc.state.spent_budget) - Float64(uc.metadata.total_budget)
        pct = (overspend / Float64(uc.metadata.total_budget)) * 100
        return (false, "spent \$$(Float64(uc.state.spent_budget)/Float64(MICRODOLLAR)) exceeds budget by $(round(pct; digits=2))%")
    end

    # Check 4: Exposure consistency
    if total_spent_exposure_rate(uc) < uc.state.spent_budget
        return (false, "totalExposure $(total_spent_exposure_rate(uc)) < spentBudget $(uc.state.spent_budget)")
    end

    # Check 5: Campaign duration valid
    if uc.metadata.start_time >= uc.metadata.end_time
        return (false, "start_time >= end_time")
    end

    # Check 6: Progress bounds
    progress = elapsed_duration_percent(uc)
    if isnan(progress) || progress < 0.0 || progress > 1.0
        return (false, "progress $progress out of range [0.0, 1.0]")
    end

    # Check 7: Pending wins reasonable (only with delayed feedback)
    if !isnothing(uc.sim_config.win_delay_generator)
        total_exposure = uc.state.spent_budget + total_budget(uc.state.pending_wins)
        if Float64(total_exposure) > Float64(uc.metadata.total_budget) * 1.10
            return (false, "total exposure (spent + pending) exceeds 110% of budget")
        end
    end

    return (true, "")
end

function check_soft_constraints(uc::SimulationUseCase)::Tuple{Bool,String}
    uc.expects.ignore_pacing_violations_while_in_flight && return (true, "")

    target_spend = round(Int64, elapsed_duration_percent(uc) * Float64(uc.metadata.total_budget))
    has_delayed_feedback = !isnothing(uc.sim_config.win_delay_generator)

    if has_delayed_feedback && uc.state.current_time < uc.metadata.end_time
        deviation = max(Int64(20) * MICRODOLLAR, round(Int64, Float64(target_spend) * 0.15))
    else
        pct = uc.expects.pacing_delta_percentage
        deviation = max(Int64(20) * MICRODOLLAR, round(Int64, Float64(target_spend) * pct))
    end

    if uc.state.spent_budget > target_spend + deviation || uc.state.spent_budget < target_spend - deviation
        return (false, "spent \$$(Float64(uc.state.spent_budget)/Float64(MICRODOLLAR)) outside target \$$(Float64(target_spend)/Float64(MICRODOLLAR)) ± \$$(Float64(deviation)/Float64(MICRODOLLAR))")
    end

    return (true, "")
end
