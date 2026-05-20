# Discrete-time MPC budget pacing simulation
# Combines adaptive Smith estimation with a safety constraint (cap at naive rate)

include("mpc_pacer.jl")
include("mpc_simulation_use_cases.jl")

# ── Main Entry Points ────────────────────────────────────────────────────────

function run_single_mpc_scenario!(uc::SimulationUseCase; plot::Bool=true, verbose::Bool=false, τ::Float64=5.0, index::Int=0)
    # Initialize
    tick_with!(uc, 0.0)

    pacer = MPCPacer(uc.state.current_time; τ)

    iterations = 0
    hard_violations = 0
    soft_violations = 0
    plot_data = PlotDataPoint[]

    while should_tick(uc, iterations)
        iterations += 1

        # Build MPC input and calculate
        result = calculate_mpc_mode!(pacer;
            current_time = uc.state.current_time,
            start_time = uc.metadata.start_time,
            end_time = uc.metadata.end_time,
            total_budget = Float64(uc.metadata.total_budget) / Float64(MICRODOLLAR),
            current_spent = Float64(total_spent_exposure_rate(uc)) / Float64(MICRODOLLAR),
            max_spend_rate = Float64(eligible_spend_budget_rate(uc)) / Float64(MICRODOLLAR) * uc.state.win_percent,
        )

        update_throttle!(uc, result.probability)

        # Hard constraint check
        passed, msg = check_hard_constraints(uc, result)
        if !passed
            hard_violations += 1
            verbose && println("  HARD VIOLATION iter $iterations: $msg")
            break
        end

        # Soft constraint check
        passed, msg = check_soft_constraints(uc)
        if !passed
            soft_violations += 1
            verbose && println("  SOFT VIOLATION iter $iterations: $msg")
        end

        # Tick
        tick_delta = tick!(uc)

        # Collect plot data
        push!(plot_data, PlotDataPoint(
            uc.state.current_time,
            tick_delta,
            linear_on_target_budget_utilization(uc),
            expected_budget_utilization(uc),
            actualized_spent_budget_percent(uc),
            uc.state.throttle,
            uc.state.last_observed_eligible_impression_rate,
        ))
    end

    verbose && println("  Completed $iterations iterations: $hard_violations hard, $soft_violations soft violations")

    # Generate plot
    if plot && !isempty(plot_data)
        plot_name = "MPCPacer_$(index)_$(uc.name)"
        plot_simulation(plot_name, plot_data, uc.metadata.start_time, uc.metadata.end_time;
            output_dir="src/mpc/plots")
    end

    return (hard_violations, soft_violations, iterations)
end

function run_mpc_simulation(; scenario::Union{Nothing,String}=nothing, plot::Bool=true, verbose::Bool=false, τ::Float64=5.0)
    use_cases = create_default_mpc_simulation_use_cases()

    if !isnothing(scenario)
        idx = findfirst(uc -> uc.name == scenario, use_cases)
        isnothing(idx) && error("Unknown scenario: $scenario. Available: $(join(map(uc -> uc.name, use_cases), ", "))")
        use_cases = [use_cases[idx]]
    end

    results = Dict{String,NamedTuple{(:hard_violations, :soft_violations, :iterations),Tuple{Int,Int,Int}}}()

    for (i, uc) in enumerate(use_cases)
        verbose && println("Running scenario $i/$(length(use_cases)): $(uc.name)")
        hv, sv, iters = run_single_mpc_scenario!(uc; plot, verbose, τ, index=i-1)
        results[uc.name] = (hard_violations=hv, soft_violations=sv, iterations=iters)
    end

    return results
end
