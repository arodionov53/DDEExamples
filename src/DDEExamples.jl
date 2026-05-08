module DDEExamples

using DelayDiffEq
using DifferentialEquations: Tsit5, ODEProblem
using Plots
using Statistics: mean, std

export solve_mackey_glass, solve_logistic_dde, solve_two_delay, solve_random_delay,
       solve_mackey_glass_nodelay, solve_logistic_nodelay, solve_two_delay_nodelay,
       solve_budget_nodelay, solve_budget_delay, solve_budget_corrected_denom, solve_budget_smith,
       solve_budget_smith_mismatch, solve_budget_smith_adaptive, demo_smith_mismatch,
       solve_budget_imc, solve_budget_mpc,
       solve_budget_pid, demo_budget_controllers, demo_budget_controllers_noise,
       solve_budget_pid_pacer, demo_pid_pacer, demo_pid_pacer_noise,
       demo_budget_delay, demo_budget_delay_with_noise,
       solve_budget_corrected_denom_spike, solve_budget_smith_spike,
       solve_budget_pid_pacer_spike, demo_demand_spike,
       solve_budget_corrected_denom_two_spikes, solve_budget_smith_two_spikes,
       solve_budget_pid_pacer_two_spikes, demo_demand_two_spikes,
       demo_demand_spike_then_drop,
       demo_demand_random,
       demo, demo_zero_delay

include("mackey_glass.jl")
include("logistic.jl")
include("two_delay.jl")
include("budget_controllers.jl")
include("imc.jl")
include("mpc.jl")
include("pid_pacer.jl")
include("demand_spike.jl")
include("demand_multi.jl")
include("demand_random.jl")
include("demos.jl")

end # module DDEExamples
