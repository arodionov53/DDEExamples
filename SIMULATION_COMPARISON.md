# Discrete-Time Simulation Comparison

Comparison of five budget pacing controllers under production-realistic stochastic conditions.

## Controllers

| Controller | Strategy | Key Idea |
|---|---|---|
| **PID Pacer** | Level-based error → PID → probability | Compares cumulative target vs total exposure; integral removes steady-state offset |
| **Smith Predictor** | Reconstruct balance via delay buffer | `B̂(t) = B(t-τ) - (S(t) - S(t-τ))`; ideal rate = predicted balance / time remaining |
| **Corrected Denominator** | Compensate observation staleness | `rate = observedBalance / (timeRemaining + τ)`; simple one-line delay compensation |
| **MPC** | Adaptive Smith + safety constraint | Estimates α̂ from observed vs predicted drain; caps rate at naive level to prevent overspend |
| **Artstein** | Integral predictor over full delay window | `B̂(t) = B(t-τ) - ∫_{t-τ}^{t} rate(s) ds`; uses trapezoidal integration of spend rate history |

## Results

All five controllers pass every scenario with **zero violations** on deterministic and moderately stochastic scenarios. The extreme scenario (`unstable_eligible_impressions_with_spikes`) can produce sporadic hard violations due to 200x traffic spikes — this is expected and the test marks it `ignore_pacing_violations_while_in_flight = true`.

### Per-Scenario Breakdown

| Scenario | PID | Smith | CorrDenom | MPC | Artstein |
|---|:---:|:---:|:---:|:---:|:---:|
| no_supply | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| exactly_1cent_per_second_spend_no_win_delay | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| exactly_1cent_per_second_spend | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_new_campaign_86400_budget | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_subsecond_intervals | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_15s_tick | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_multi_minute_intervals | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| hard_catch_up | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| resumed_near_end_of_flight | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| resumed_over_paced_near_end_of_flight | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| compensate_for_spikes | 0H 0S | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| unstable_eligible_impressions_with_spikes | 0–1H* | 0–1H* | 0H 0S | 0–1H* | 0–1H* |

\* Stochastic — 200x traffic spikes with low base rate (10 imp/s, σ=200) can trigger budget overshoot regardless of controller. This scenario has `ignore_pacing_violations_while_in_flight = true`.

### Unstable Scenario Reliability (5 runs each)

| Controller | Hard violations across 5 runs |
|---|---|
| PID | [0, 0, 0, 0, 0] |
| Smith | [0, 1, 0, 0, 0] |
| CorrDenom | [0, 0, 0, 0, 0] |
| MPC | [0, 0, 1, 0, 0] |
| Artstein | [1, 1, 0, 0, 1] |

CorrDenom and PID appear most robust to extreme spikes, followed by Smith and MPC. Artstein is slightly more prone to overshoot under bursty conditions because the integral predictor accumulates noisy rate estimates.

## Constraint Definitions

**Hard constraints** (violation = immediate abort):
- Throttle must be 0 outside campaign flight window
- Throttle must be in [0.0, 1.0]
- Spent budget must not exceed 105% of total budget
- Total exposure ≥ spent budget (consistency)
- Campaign duration valid (start < end)
- Progress in [0.0, 1.0]
- Total exposure (spent + pending) must not exceed 110% of budget

**Soft constraints** (violation = tracked but not fatal):
- Spent budget must track linear target within ±15% (with delayed feedback) or ±configured delta (without)

## Simulation Environment

Each scenario runs under production-realistic conditions:
- **Stochastic supply**: eligible impressions drawn from Normal(μ, σ) with occasional 200x spikes
- **Auction noise**: win rate ~25% ± 3% (Normal), varying per tick
- **Delayed win notifications**: bimodal latency distribution (~2 min modal peak), calibrated from 11.5B real events
- **Stale supply observations**: eligible impression rate arrives 30–300s late (log-normal)
- **Variable tick intervals**: sub-second to multi-minute, drawn from Normal(μ, σ)

## Scenarios

| # | Scenario | What It Tests |
|---|---|---|
| 1 | no_supply | Zero eligible impressions — controller must output 0 |
| 2 | exactly_1cent_per_second_spend_no_win_delay | Perfect conditions, no feedback delay |
| 3 | exactly_1cent_per_second_spend | Perfect supply but with delayed win notifications |
| 4 | pace_new_campaign_86400_budget | Standard campaign: $86,400 over 1 hour |
| 5 | pace_with_fuzzy_subsecond_intervals | Very frequent ticks (0.5s ± 0.5s) |
| 6 | pace_with_fuzzy_15s_tick | Moderate tick rate (15s ± 1s) |
| 7 | pace_with_fuzzy_multi_minute_intervals | Slow ticks (15s, no jitter) |
| 8 | hard_catch_up | 25-hour campaign, simulation starts in last hour — must catch up |
| 9 | resumed_near_end_of_flight | Pre-loaded 90% spent, must pace remaining 10% |
| 10 | resumed_over_paced_near_end_of_flight | Pre-loaded 95% spent, must throttle down |
| 11 | compensate_for_spikes | 3-day campaign with 0.2% chance of 200x traffic spikes |
| 12 | unstable_eligible_impressions_with_spikes | Extreme: low base rate (10), high stddev (200), frequent spikes |

## Controller Design Comparison

| Aspect | PID | Smith | CorrDenom | MPC | Artstein |
|---|---|---|---|---|---|
| **Delay compensation** | Implicit (integral accumulates) | Point-difference S(t) - S(t-τ) | Adds τ to denominator | Adaptive α̂ + Smith | ∫ rate(s) ds over [t-τ, t] |
| **History used** | Last error only | 1 point at t-τ | None (stateless) | 2 points at t-τ, t-2τ | Full buffer over [t-τ, t] |
| **Model complexity** | 3 tunable params (Kp, Ki, Kd) | 1 param (τ) | 1 param (τ) | 1 param (τ) + adaptive α̂ | 1 param (τ) |
| **Robustness to noise** | High (integral smooths) | Moderate | High (simple formula) | Moderate (α̂ can be noisy) | Moderate (integral of noisy rate) |
| **Theoretical optimality** | Sub-optimal (5–20% under-delivery in DDE model) | Optimal when model is exact | Near-optimal at small τ | Optimal + safe when α̂ accurate | Optimal (equivalent to Smith in continuous limit) |
| **Safety under uncertainty** | Conservative (never overspends) | Can overspend if model wrong | Conservative (under-delivers at large τ) | Explicitly capped at safe rate | Can overspend under extreme noise |

## Analysis

In the discrete-time simulation all five controllers are equally **safe** under normal-to-moderate stochastic conditions (scenarios 1–11). The `totalExposure` signal (spent + pending + projected outstanding) provides partial delay compensation that benefits all controllers equally.

Under **extreme conditions** (scenario 12: 200x spikes with near-zero baseline), the controllers rank by robustness:
1. **CorrDenom** — most robust; simple formula has no state that can accumulate errors
2. **PID** — robust; integral provides damping against spike-induced oscillations
3. **Smith/MPC** — occasional violations when spike timing aligns poorly with the delay window
4. **Artstein** — slightly more fragile; integrating noisy rates over the full window can amplify spike effects

The theoretical differences (Smith = optimal, PID = conservative) are more visible in the continuous-time DDE models. See `dde-vs-go-pid-simulation.md` for the full comparison.

## How to Reproduce

```bash
./refresh-simulation-charts.sh
```

Or individually:

```julia
using DDEExamples
run_pid_simulation(; verbose=true)
run_smith_simulation(; verbose=true)
run_corrdenom_simulation(; verbose=true)
run_mpc_simulation(; verbose=true)
run_artstein_simulation(; verbose=true)
```

Run tests:

```bash
julia --project=. test/runtests.jl
```

Plots are saved to:
- `src/pid/plots/PIDPacer_*`
- `src/smith/plots/SmithPacer_*`
- `src/corrdenom/plots/CorrDenomPacer_*`
- `src/mpc/plots/MPCPacer_*`
- `src/artstein/plots/ArtsteinPacer_*`
