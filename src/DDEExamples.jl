module DDEExamples

using DelayDiffEq
using DifferentialEquations: Tsit5, ODEProblem
using Plots
using Statistics: mean, std

# ── Shared simulation infrastructure ─────────────────────────────────────────

const MICRODOLLAR = Int64(1_000_000)
const MICROCENT = MICRODOLLAR ÷ 100
const MICROMILLI = Int64(1000)

include("simulation/pending_queues.jl")
include("simulation/simulation.jl")
include("simulation/generators.jl")
include("simulation/plot_simulation.jl")

# ── Controller-specific simulations ─────────────────────────────────────────

include("pid/pid_simulation.jl")
include("smith/smith_simulation.jl")

export solve_mackey_glass, solve_logistic_dde, solve_two_delay, solve_random_delay,
       solve_mackey_glass_nodelay, solve_logistic_nodelay, solve_two_delay_nodelay,
       solve_budget_nodelay, solve_budget_delay, solve_budget_corrected_denom, solve_budget_smith,
       solve_budget_smith_mismatch, solve_budget_smith_adaptive, demo_smith_mismatch,
       solve_budget_imc, solve_budget_mpc,
       solve_budget_mpc_spike, solve_budget_imc_spike, demo_mpc_spike,
       solve_budget_mpc_two_spikes, solve_budget_imc_two_spikes, demo_mpc_two_spikes,
       solve_budget_pid, demo_budget_controllers, demo_budget_controllers_noise,
       solve_budget_pid_pacer, demo_pid_pacer, demo_pid_pacer_noise,
       demo_budget_delay, demo_budget_delay_with_noise,
       solve_budget_corrected_denom_spike, solve_budget_smith_spike,
       solve_budget_pid_pacer_spike, demo_demand_spike,
       solve_budget_corrected_denom_two_spikes, solve_budget_smith_two_spikes,
       solve_budget_pid_pacer_two_spikes, demo_demand_two_spikes,
       demo_demand_spike_then_drop,
       demo_demand_random,
       demo, demo_zero_delay,
       run_pid_simulation, run_single_scenario!,
       create_default_pid_pacer_simulation_use_cases,
       PIDPacer, PIDConfig, PIDResult, calculate_cruise_mode!,
       run_smith_simulation, run_single_smith_scenario!,
       create_default_smith_simulation_use_cases,
       SmithPacer, SmithConfig, SmithResult, calculate_smith_mode!,
       SimulationUseCase, CampaignMetadata, CampaignState, SimulationConfig,
       tick!, tick_with!, should_tick,
       plot_simulation

include("mackey_glass.jl")
include("logistic.jl")
include("two_delay.jl")
include("budget_controllers.jl")
include("imc.jl")
include("mpc.jl")
include("mpc_spike.jl")
include("pid_pacer.jl")
include("demand_spike.jl")
include("demand_multi.jl")
include("demand_random.jl")
include("demos.jl")

end # module DDEExamples
