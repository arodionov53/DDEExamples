# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Activate project and install deps
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run tests (discrete-time PID simulation scenarios)
julia --project=. test/runtests.jl

# Regenerate all DDE model comparison charts (continuous-time)
./refresh-ddeexamples-charts.sh

# Regenerate discrete-time simulation plots (PID + Smith)
./refresh-simulation-charts.sh

# Run a single simulation interactively
julia --project=. -e 'using DDEExamples; run_pid_simulation(; scenario="pace_new_campaign_86400_budget", verbose=true)'
julia --project=. -e 'using DDEExamples; run_smith_simulation(; scenario="pace_new_campaign_86400_budget", verbose=true)'
```

## Architecture

This package explores budget pacing under information delay from two complementary angles:

### 1. Continuous-time DDE models (`src/DDEExamples.jl`)

The main module file contains all DDE solver functions and demo/plotting functions. These model budget controllers as continuous-time Delay Differential Equations solved with `MethodOfSteps(Tsit5())`. Four controllers are compared: Naive, Corrected denominator, Smith predictor, and PID pacer.

Demand-spike variants use a segment-stitching pattern: solve up to the discontinuity, then restart the solver with modified initial conditions and a history function that replays prior segments.

### 2. Discrete-time simulation (`src/simulation/` + `src/pid/` + `src/smith/`)

A faithful port of the Go production PID simulation (`pacing-oracle/pkg/pacing/simulation/`). This tests the actual controller logic under stochastic market conditions.

**Simulation engine** (`src/simulation/`):
- `simulation.jl` — `SimulationUseCase` struct hierarchy, `tick!`/`tick_with!` engine, helper methods
- `pending_queues.jl` — delayed win notification and rate observation delivery queues
- `generators.jl` — stochastic generators for supply, win rate, tick intervals, and latency distributions (calibrated from 11.5B real events)
- `plot_simulation.jl` — per-scenario plot generation
- `check_constraints.jl` — hard/soft invariant checks

**Controllers**:
- `src/pid/pid_pacer.jl` — `PIDPacer` struct and `calculate_cruise_mode!` (level-based error: targetSpent - totalExposure)
- `src/smith/smith_pacer.jl` — `SmithPacer` struct and `calculate_smith_mode!` (reconstructs balance via delay buffer)
- `src/pid/pid_simulation.jl` + `src/smith/smith_simulation.jl` — scenario runners (`run_single_scenario!`, `run_single_smith_scenario!`)
- `src/pid/pid_simulation_use_cases.jl` + `src/smith/smith_simulation_use_cases.jl` — test scenario definitions

### Key design decisions

- All monetary values use **microdollars** (`Int64`, 1 dollar = 1_000_000 microdollars) to avoid floating-point error. Constants: `MICRODOLLAR`, `MICROCENT`, `MICROMILLI`.
- Plot naming convention: `PacerType_Index_ScenarioName` (e.g., `PIDPacer_3_pace_new_campaign_86400_budget`).
- The Go PID uses level-based error (cumulative $) while the Julia DDE model uses rate-based error ($/second). See `dde-vs-go-pid-simulation.md` for the full comparison.
- Smith pacer normalization: `probability = idealRate / maxSpendRate` (no dt multiplication), matching the Go implementation.

## Reference implementation

The Julia discrete-time simulation is aligned with the Go implementation at `pacing-oracle/pkg/pacing/simulation/`. When modifying controller logic or simulation mechanics, the Go code is the source of truth.
