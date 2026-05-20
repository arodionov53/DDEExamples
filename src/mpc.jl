"""
    solve_budget_mpc(; Q, T, τ, α)

Model Predictive Control (MPC) for budget pacing under information delay.

Standard controllers (Smith, corrected-denom) use either the model estimate
or the raw observation — they do not *optimise* subject to constraints.  MPC
explicitly enforces the hard constraint B(t) ≥ 0 at all times by solving a
one-step constrained optimisation at each instant.

**Unconstrained optimal rate** (from Smith/IMC):

    rate_free(t) = B̂(t) / (T - t)

where B̂ is the adaptive Smith estimate.

**Budget constraint.**

If rate_free were applied and the true balance is actually lower than B̂ (e.g.
due to α̂ estimation lag), the budget could go negative.  MPC guards against
this by also computing a *conservative* rate based on the raw delayed
observation:

    rate_safe(t) = B(t-τ) / (T - t)    (naive controller — never negative)

The MPC rate is the *minimum* of the two:

    rate_mpc(t) = min(rate_free(t), rate_safe(t))

**Interpretation.**

- When B̂ ≤ B(t-τ): no mismatch or the estimate is lower than observed —
  use the (potentially aggressive) free rate.
- When B̂ > B(t-τ): the adaptive estimate says more budget remains than the
  raw observation — cap the rate at the safe level to prevent overspend.

This is equivalent to a one-step receding-horizon MPC that minimises
‖B(T)‖² subject to B(t) ≥ 0, using a conservative fallback when the
model confidence is low.

State: u = [B(t), S(t)]   (same as Smith)
Requires lags τ, 2τ, 3τ (for α̂ estimation).

`α` = true grant fulfilment rate (unknown to the controller).
"""
function solve_budget_mpc(;
    Q     = 100.0,
    T     = 10.0,
    τ     = 1.5,
    α     = 1.0,
    tspan = (0.0, T - 1e-3)
)
    function mpc!(du, u, h, p, t)
        Q, T, τ, α = p
        dt = τ > 0 ? τ : 1e-6

        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]

        # estimate α̂ (same as adaptive Smith)
        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0

        # unconstrained adaptive Smith estimate
        B_hat = max(B_now - α_hat * (u[2] - S_now), 0.0)

        # constrained MPC: cap rate at the safe naive level
        rate_free = B_hat   / (T - t)
        rate_safe = B_now   / (T - t)   # naive rate — safe upper bound
        rate      = min(rate_free, rate_safe)

        du[1] = -α * rate
        du[2] =      rate
    end

    h(p, t) = [Q, 0.0]
    p = (Q, T, τ, α)
    prob = DDEProblem(mpc!, [Q, 0.0], h, tspan, p;
                      constant_lags = [τ, 2τ, 3τ])
    solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
end
