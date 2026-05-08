# DDEExamples

A Julia package with worked examples of **Delay Differential Equations (DDEs)** solved with [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/).

Each example pairs a real-world motivation with a runnable solver function and a comparison against the zero-delay (ODE) limit.  Detailed results and explanations are in:

- [RESULTS.md](RESULTS.md) — classical DDEs (Mackey-Glass, logistic, two-delay, zero-delay comparison)
- [RESULTS_BUDGET.md](RESULTS_BUDGET.md) — budget-pacing controllers (sections 5–5g: naive, corrected denom, Smith, PID, PIDPacer, demand spikes, random events)

## Examples

| # | Name | Equation | Key behavior |
|---|------|----------|-------------|
| 1 | Mackey-Glass | `du/dt = β u(t-τ)/(1+u(t-τ)ⁿ) - γ u` | Deterministic chaos |
| 2 | Delayed logistic growth | `du/dt = r u(t) (1 - u(t-τ)/K)` | Boom-and-bust oscillations |
| 3 | Two-delay system | `du/dt = -a u(t-τ₁) - b u(t-τ₂)` | Damped oscillations |
| 4 | Random delay logistic | Same as #2, τ ~ Uniform[τ_min, τ_max] | Monte Carlo ensemble |
| 5 | Budget spending with delay | `dB/dt = -B(t-τ)/(T-t)` | Overspend from stale information |

Example 5 includes four controllers of increasing sophistication:

- **Naive** — `dB/dt = -B(t-τ)/(T-t)`, baseline; systematically under-spends early and rushes near the deadline
- **Corrected denominator** — `dB/dt = -B(t-τ)/(T-t+τ)`, a one-line fix that works well for small delays
- **Smith predictor** — reconstructs the true current balance from the delayed observation and cumulative spend history; achieves near-perfect pacing at any delay
- **PID** — tracks the reference trajectory `B_ref(t) = Q(1-t/T)` via a delayed error signal with integral and derivative terms

### Production PIDPacer (section 5c)

`solve_budget_pid_pacer` / `demo_pid_pacer` model the production Go implementation faithfully:
- Rate-based error signal: `e(t) = targetRate − observedRate(t−τ)`
- Observed rate estimated from finite differences of delayed balance
- PID output mapped through a sigmoid to a grant probability
- Integral windup protection via clamping

### Noise sensitivity (sections 5b, 5c)

`demo_budget_delay_with_noise` and `demo_pid_pacer_noise` compare all strategies across multiple τ noise levels (mean ± 1σ ribbons):

![Effect of τ noise on spending strategies](plots/budget_delay_noise.png)

| Strategy | τ = 10%·T | τ = 50%·T | Noise sensitivity |
|----------|-----------|-----------|-------------------|
| Naive | 111% spent | 413% spent | Low bias change, moderate spread |
| Corrected denom | 100% | 95% (under-spends) | Bias dominates, noise negligible |
| Smith predictor | 100% | 100% | Immune — all curves overlap |

### Demand spikes (sections 5d–5g)

Each controller is also stress-tested against instantaneous demand events — budget discontinuities with duration < τ, so the controller is blind to them for a full delay window.

| Demo function | Scenario |
|---------------|----------|
| `demo_demand_spike` | Single spike at `t_spike_frac·T` |
| `demo_demand_two_spikes` | Two spikes; second may arrive before recovery from first |
| `demo_demand_spike_then_drop` | Spike followed by a budget refund; tests windfall exploitation |
| `demo_demand_random` | Monte Carlo: N random signed events, mean ± 1σ ribbon |

## Installation

```julia
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
```

## Quick start

```julia
using DDEExamples

# Run all classic DDE examples and save dde_examples.png
DDEExamples.demo()

# Compare DDEs with their zero-delay (ODE) counterparts
DDEExamples.demo_zero_delay()

# Budget spending strategies under information delay
demo_budget_delay()                          # default: 5%, 10%, 50% of T
demo_budget_delay(delays = [1.0, 2.0, 3.0]) # custom delays
demo_budget_delay(tau_noise = 0.03)          # add ±3% noise to τ

# All four controllers side-by-side (saves budget_controllers.png)
demo_budget_controllers()
demo_budget_controllers_noise()              # with τ uncertainty bands

# Noise sensitivity for three core strategies (saves budget_delay_noise.png)
demo_budget_delay_with_noise()
demo_budget_delay_with_noise(noise_levels = [0.0, 0.1, 0.5])

# Production PIDPacer model (saves pid_pacer_comparison.png / pid_pacer_noise.png)
demo_pid_pacer()
demo_pid_pacer_noise()

# Demand-spike experiments
demo_demand_spike()                          # single spike, duration < τ
demo_demand_two_spikes()                     # two successive spikes
demo_demand_spike_then_drop()                # spike then budget refund
demo_demand_random()                         # Monte Carlo random events
```

### Solve individually

```julia
sol = solve_mackey_glass(τ = 4.0)        # chaotic Mackey-Glass
sol = solve_logistic_dde(τ = 5.0)        # boom-and-bust logistic
sol = solve_two_delay(τ₁ = 1.0, τ₂ = 3.0)
sim = solve_random_delay(trajectories = 50)  # ensemble solution

sol = solve_budget_delay(Q = 100.0, T = 10.0, τ = 1.5)          # naive
sol = solve_budget_corrected_denom(Q = 100.0, T = 10.0, τ = 1.5) # corrected denom
sol = solve_budget_smith(Q = 100.0, T = 10.0, τ = 1.5)           # Smith predictor
sol = solve_budget_pid(Q = 100.0, T = 10.0, τ = 1.5)             # PID
sol = solve_budget_pid_pacer(Q = 100.0, T = 10.0, τ = 1.0)       # production PIDPacer
```

All solvers return a `DifferentialEquations.jl` solution object:

```julia
sol.t          # time points
sol[1, :]      # state variable values
sol(3.7)       # interpolate at any time
```

## Solver

All examples use `MethodOfSteps(Tsit5())` — the standard Method of Steps wrapping a 5th-order adaptive Runge-Kutta solver.

## References

- [DDEProblem API](https://docs.sciml.ai/DiffEqDocs/stable/types/dde_types/) — constructor, history function interface, `constant_lags` vs `dependent_lags`

## Dependencies

- [DelayDiffEq.jl](https://github.com/SciML/DelayDiffEq.jl)
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/) (stdlib)
