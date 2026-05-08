# Discrete-time PID budget pacing simulation
# Faithful port of the Go simulation at pacing-oracle/pkg/pacing/simulation/

# ── Constants ────────────────────────────────────────────────────────────────

const MICRODOLLAR = Int64(1_000_000)
const MICROCENT = MICRODOLLAR ÷ 100
const MICROMILLI = Int64(1000)

# ── Core Types ───────────────────────────────────────────────────────────────

include("pid_pacer.jl")
include("../simulation/pending_queues.jl")
include("../simulation/simulation.jl")

# ── Stochastic Generators ────────────────────────────────────────────────────

include("../simulation/generators.jl")

# ── Test Scenarios ───────────────────────────────────────────────────────────

include("pid_simulation_use_cases.jl")

# ── Main Entry Points ────────────────────────────────────────────────────────

function run_single_scenario!(uc::SimulationUseCase; plot::Bool=true, verbose::Bool=false)
    # Initialize
    tick_with!(uc, 0.0)

    pacer = PIDPacer(uc.state.current_time)

    iterations = 0
    hard_violations = 0
    soft_violations = 0
    plot_data = PlotDataPoint[]

    while should_tick(uc, iterations)
        iterations += 1

        # Target spent: linear over campaign duration
        target_spent = Float64(uc.metadata.total_budget) * elapsed_duration_percent(uc)

        # Build PID input and calculate
        result = calculate_cruise_mode!(pacer;
            current_time = uc.state.current_time,
            start_time = uc.metadata.start_time,
            end_time = uc.metadata.end_time,
            target_spent = target_spent / Float64(MICRODOLLAR),
            total_exposure = Float64(total_spent_exposure_rate(uc)) / Float64(MICRODOLLAR),
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
        plot_simulation(uc.name, plot_data, uc.metadata.start_time, uc.metadata.end_time)
    end

    return (hard_violations, soft_violations, iterations)
end

function run_pid_simulation(; scenario::Union{Nothing,String}=nothing, plot::Bool=true, verbose::Bool=false)
    use_cases = create_default_pid_pacer_simulation_use_cases()

    if !isnothing(scenario)
        idx = findfirst(uc -> uc.name == scenario, use_cases)
        isnothing(idx) && error("Unknown scenario: $scenario. Available: $(join(map(uc -> uc.name, use_cases), ", "))")
        use_cases = [use_cases[idx]]
    end

    results = Dict{String,NamedTuple{(:hard_violations, :soft_violations, :iterations),Tuple{Int,Int,Int}}}()

    for (i, uc) in enumerate(use_cases)
        verbose && println("Running scenario $i/$(length(use_cases)): $(uc.name)")
        hv, sv, iters = run_single_scenario!(uc; plot, verbose)
        results[uc.name] = (hard_violations=hv, soft_violations=sv, iterations=iters)
    end

    return results
end

# ── Plotting ─────────────────────────────────────────────────────────────────

include("../simulation/plot_simulation.jl")
