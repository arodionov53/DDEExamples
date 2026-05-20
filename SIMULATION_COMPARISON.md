# Discrete-Time Simulation Comparison

Comparison of four budget pacing controllers under production-realistic stochastic conditions.

## Controllers

| Controller | Strategy | Key Idea |
|---|---|---|
| **PID Pacer** | Level-based error → PID → probability | Compares cumulative target vs total exposure; integral removes steady-state offset |
| **Smith Predictor** | Reconstruct balance via delay buffer | `B̂(t) = B(t-τ) - (S(t) - S(t-τ))`; ideal rate = predicted balance / time remaining |
| **Corrected Denominator** | Compensate observation staleness | `rate = observedBalance / (timeRemaining + τ)`; simple one-line delay compensation |
| **MPC** | Adaptive Smith + safety constraint | Estimates α̂ from observed vs predicted drain; caps rate at naive level to prevent overspend |

## Results

All four controllers pass every scenario with **zero hard and zero soft violations**.

### Per-Scenario Breakdown

| Scenario | PID | Smith | CorrDenom | MPC |
|---|:---:|:---:|:---:|:---:|
| no_supply | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| exactly_1cent_per_second_spend_no_win_delay | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| exactly_1cent_per_second_spend | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_new_campaign_86400_budget | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_subsecond_intervals | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_15s_tick | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| pace_with_fuzzy_multi_minute_intervals | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| hard_catch_up | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| resumed_near_end_of_flight | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| resumed_over_paced_near_end_of_flight | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| compensate_for_spikes | 0H 0S | 0H 0S | 0H 0S | 0H 0S |
| unstable_eligible_impressions_with_spikes | 0H 0S | 0H 0S | 0H 0S | 0H 0S |

### Constraint Definitions

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

## Analysis

All four controllers are equally **safe** in the discrete-time simulation — none violates budget or pacing constraints under any tested scenario. This is because the simulation's `totalExposure` signal (spent + pending + projected outstanding) provides partial delay compensation that all controllers benefit from.

The theoretical differences between controllers are more visible in the **continuous-time DDE models** (see `dde-vs-go-pid-simulation.md`), where the controller sees only a purely delayed balance signal:

| Controller | Continuous-Time Behavior |
|---|---|
| **Smith** | Perfect pacing at any delay (mathematically optimal) |
| **MPC** | Same as Smith when α̂ is accurate; safely conservative when uncertain |
| **Corrected Denominator** | Near-perfect at small delays; slight under-delivery at large delays |
| **PID** | Never overspends but systematically under-delivers by 5–20% |

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
```

Plots are saved to:
- `src/pid/plots/PIDPacer_*`
- `src/smith/plots/SmithPacer_*`
- `src/corrdenom/plots/CorrDenomPacer_*`
- `src/mpc/plots/MPCPacer_*`
