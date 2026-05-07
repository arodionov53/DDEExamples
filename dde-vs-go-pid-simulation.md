# DDE (Julia) vs Go PID Simulation

## What This Simulation System Does

You have an ad campaign budget pacing problem: spend $X evenly over a time window, but you're making decisions based on *stale* data. You're exploring that problem from two angles.

---

## Julia (DDEExamples) — "What does the math say?"

The Julia code models the pacing problem as a **continuous-time Delay Differential Equation**. It strips away all the production messiness (stochastic supply, auctions, batching) and asks the pure control-theory question:

> If a controller can only see its balance from τ seconds ago, what spending trajectory does it produce?

It compares four control strategies:

| Controller | How it works | Result |
|---|---|---|
| **Naive** | Spend based on stale balance ÷ remaining time | Catastrophically overspends at large delays |
| **Corrected denominator** | Add τ to the denominator to pre-compensate for staleness | Near-perfect at small delays, slight error at large ones |
| **Smith predictor** | Reconstruct true balance using `B(t-τ) - (spent_since_t-τ)` | Always perfect, any delay |
| **PID Pacer** (models the Go code) | Rate-based error → PID → probability | Never overspends but systematically under-delivers by 5-20% |

It then subjects all four to: random delay noise, demand spikes (sudden external spend), budget refunds, and sequences of random events. The Smith predictor is immune to all of them. The PID pacer is *safe* (never blows the budget) but *conservative* (leaves money on the table).

**Purpose:** Identify the theoretical performance ceiling and failure modes of the chosen PID approach, and motivate whether a Smith predictor component is worth adding.

---

## Go (simulation package) — "Does it survive production reality?"

The Go code runs the **actual production PID controller** (`pid_pacer.go`) inside a discrete-time simulation that models everything the Julia code deliberately ignores:

1. **Random supply** — eligible impressions/second drawn from a normal distribution, with occasional 200x spikes
2. **Auction wins** — you bid but only win ~25% of auctions (varying randomly)
3. **Delayed win notifications** — modeled with empirically-calibrated latency distributions measured from 11.5B real events:
   - Load events: ~0.5s mean (exponential)
   - Impression events: ~2 minute modal peak (bimodal mixture)
   - Tracker events: ~4 minute median (log-normal, heavy tail)
   - Click events: multi-modal across 20s, 2m, 1h, 12h (unusable for closed-loop control)
4. **Stale supply observations** — the eligible impression rate the controller sees is 30-300s old
5. **Variable tick intervals** — the controller fires at jittery intervals (sub-second to multi-minute)

### The loop each tick

1. Advance clock by a random interval
2. Deliver any queued win notifications whose delay has elapsed
3. Deliver any queued supply-rate observations whose delay has elapsed
4. Generate new ground-truth supply and win-rate from random distributions
5. Compute actual wins = throttle × true_supply × win_rate × tick_duration
6. Queue those wins for delayed delivery (split into batches with random delays)
7. Feed the *stale observed state* into the real PID controller
8. PID outputs a new throttle (probability 0-1)
9. Check hard invariants (never >105% budget, throttle in bounds, progress valid) and soft invariants (spending tracks linear target within 15%)

### Test cases

The test cases cover: zero supply, exactly-matching supply, high-budget campaigns, sub-second ticks, 15-second ticks, multi-minute ticks, catch-up after late start, resume when already 90% or 95% spent, traffic spikes, and unstable/bursty supply patterns. Each generates a gnuplot PNG for visual inspection.

**Purpose:** Prove the production controller doesn't violate safety constraints under realistic stochastic conditions, and generate visual evidence of how it behaves.

---

## How They Fit Together

| Dimension | Julia | Go |
|---|---|---|
| Time model | Continuous (adaptive DDE solver) | Discrete ticks |
| Supply | Constant or absent | Stochastic, spiky |
| Delay | Fixed τ (parameter) | Empirical distribution per event |
| Controller tested | Mathematical model of the Go PID | The actual Go PID code |
| What it proves | Theoretical limits and comparison to optimal | No hard constraint violations in practice |
| Output | Comparison plots showing trajectory vs ideal | Per-scenario gnuplot + pass/fail invariant checks |

The Julia work identified that the PID design is fundamentally conservative (under-delivers) and that a Smith predictor would be optimal. The Go simulation validates that "conservative" translates to "safe" under real-world chaos — you won't overspend, even with 200x traffic spikes, 2-minute feedback delays, and variable tick rates.

---

## Control Loop Comparison

Both systems implement the same PID control logic, but they differ in how the feedback loop is structured, what information the controller sees, and how the plant (the thing being controlled) responds.

### Error Signal

| Aspect | Julia DDE | Go Simulation |
|---|---|---|
| **Error definition** | `e = targetRate - observedRate(t-τ)` | `e = targetSpent - totalExposure` |
| **What it means** | "Am I spending at the right *rate*?" | "Have I spent the right *amount* so far?" |
| **Units** | $/second (flow) | $ (cumulative stock) |

The Julia model uses a **rate-based** error: it estimates the current spend rate from a finite difference of delayed balance observations (`(B(t-2τ) - B(t-τ)) / τ`) and compares that to the target rate `Q/T`. This makes the delay explicit — the controller literally cannot see anything newer than τ seconds ago.

The Go production code uses a **level-based** error: it compares where cumulative spending *should* be (linear target) against where it *appears* to be (total exposure = realized spend + pending wins + projected outstanding). The delay is implicit — it hides inside the staleness of `SpentBudget` (which only reflects delivered win notifications) and the `PendingWins` queue (which estimates in-flight spend).

### What the Controller Sees

| Input | Julia DDE | Go Simulation |
|---|---|---|
| Current balance | `B(t-τ)` — exactly τ seconds stale | `SpentBudget` — stale by however long win notifications take to arrive (empirical distribution, typically 30s-4min) |
| Supply rate | Not modeled (constant `request_rate = Q/T`) | `LastObservedEligibleImpressionRate` — stale by 30-300s due to `EligibleImpressionLatencyGenerator` |
| Win rate | Not modeled (implicit in `request_rate`) | `WinPercent` — observed directly each tick (no delay on this signal) |
| Pending exposure | Not modeled | `PendingWins.TotalBudget()` — the controller knows about in-flight spend it hasn't yet received confirmation for |

The Julia model is deliberately minimal: one delay parameter τ controls everything. The Go simulation layers multiple independent delays (win notifications, supply observations) on top of stochastic variation in the signals themselves.

### PID Calculation

**Julia:**
```
observed_rate = (B(t-2τ) - B(t-τ)) / τ
e             = target_rate - observed_rate
I_new         = clamp(I + e, -max_integral, max_integral)
de            = (e - e_prev) / τ
pid_out       = Kp·e + Ki·I_new + Kd·de
probability   = clamp(pid_out / (request_rate · τ), 0, 1)
spend_rate    = probability × request_rate
```

**Go:**
```
e             = targetSpent - totalExposure
integral     += e × dt
integral      = clamp(integral, -max_integral, max_integral)
derivative    = (e - lastError) / dt
CV            = Kp·e + Ki·integral + Kd·derivative
probability   = clamp(CV / (maxSpendRate × dt), 0, 1)
```

Key differences:

1. **Integral accumulation**: Julia integrates the raw error `e` (dimensionless rate error). Go integrates `e × dt` (dollar-seconds), making the integral scale with how long the error persists.

2. **Derivative source**: Julia computes derivative from two successive rate errors separated by τ. Go computes derivative from the change in cumulative error between consecutive ticks (separated by variable `dt`).

3. **Output normalization**: Both divide PID output by `maxSpendRate × dt` to convert from "dollars of correction needed" to "probability of granting the next request." Julia uses the fixed `request_rate × τ`; Go uses the *observed* `maxSpendRate × dt`, which varies with supply and win rate.

4. **dt vs τ**: In Julia, the time constant is always τ (the fixed delay). In Go, `dt` is the actual elapsed wall-clock time since the last calculation — it varies per tick and is independent of the feedback delay.

### Plant Response (How Spending Actually Happens)

**Julia:**
```
dB/dt = -probability × request_rate
```

The "plant" is trivial: spending happens continuously and deterministically at the rate the controller dictates. There is no randomness, no auction, no batching. The only complication is that the controller sees `B(t-τ)` instead of `B(t)`.

**Go:**
```
actual_impressions = throttle × actual_eligible_rate × tick_duration × win_percent
actual_budget      = actual_impressions × cost_per_impression
→ queued for delivery after random delay (ImpressionEventLatencyGenerator)
```

The plant is stochastic and multi-stage:
- The controller sets a throttle (probability)
- The market provides random supply (`actual_eligible_rate ≠ observed_rate`)
- Auctions produce random wins (`win_percent` varies)
- Win notifications arrive after random delays (log-normal, 2-minute modal peak)
- Only *delivered* notifications update `SpentBudget`

This means the Go controller is fighting three sources of uncertainty simultaneously: supply noise, auction noise, and feedback delay noise. The Julia controller only fights one: a single fixed delay τ.

### Feedback Topology

```
Julia (single-loop, single-delay):

  target_rate ──→ [error] ──→ [PID] ──→ probability ──→ [plant: dB/dt = -p·r]
                    ↑                                           │
                    └──── observed_rate ← finite_diff(B(t-τ)) ←─┘


Go (multi-signal, multi-delay):

  target_spent ──→ [error] ──→ [PID] ──→ throttle ──→ [market: supply × win% × throttle]
                     ↑                                         │
                     │                              ┌──────────┘
                     │                              ↓
                     │                    [pending win queue] ──delay──→ SpentBudget
                     │                              │
                     ├── totalExposure ←── SpentBudget + PendingWins + ProjectedOutstanding
                     │
                     └── maxSpendRate ←── ObservedEligibleRate (delayed 30-300s) × WinPercent
```

The Go system has a richer feedback structure:
- **TotalExposure** combines three signals at different staleness levels: realized spend (most stale), pending wins (partially stale), and projected outstanding (estimate based on current throttle × observed rate).
- **MaxSpendRate** is itself a delayed signal, creating a second feedback path where supply staleness affects the PID normalization denominator.
- The Julia model collapses all of this into a single delay parameter τ applied to one signal.

### Implications

The Julia model reveals the *fundamental* limitation: any PID controller operating on a delayed signal will under-deliver because the integral takes time to ramp up, and the delay prevents immediate correction. This is inherent to the control architecture regardless of implementation details.

The Go model reveals the *practical* behavior: the multi-signal feedback structure (exposure = spent + pending + projected) partially compensates for the delay by giving the controller a forward-looking estimate. The `ProjectedOutstandingSpentBudgetRate` term (`throttle × observedRate × costPerImpression`) acts as a crude predictor — not as mathematically rigorous as a Smith predictor, but it provides some delay compensation that the Julia model's pure `B(t-τ)` observation does not have.

This explains why the Go simulation often performs *better* than the Julia DDE model predicts at equivalent delays: the production controller's TotalExposure signal is less stale than a pure τ-delayed balance observation because it includes the pending-win estimate and the projected-outstanding term.
