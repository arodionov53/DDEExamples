using DDEExamples

# Run all three examples and save a combined plot
DDEExamples.demo()

# Or solve individually and inspect:
#
#   sol = solve_mackey_glass(τ = 4.0)    # increase delay → more chaotic
#   sol = solve_logistic_dde(r = 1.5)    # faster growth rate → oscillations
#   sol = solve_two_delay(τ₁ = 2.0)      # shift first delay
#
# Each function returns a DifferentialEquations.jl solution object, so you
# can index it (sol[1, :]), interpolate it (sol(15.3)), or plot it directly.
