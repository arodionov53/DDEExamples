"""
    solve_budget_imc(; Q, T, τ, λ, α)

Internal Model Control (IMC) with robustness filter for budget pacing.

**How IMC differs from Smith predictor.**

The Smith predictor computes an exact balance estimate using the grant tally:

    B̂_smith = B(t-τ) - α̂ · (S(t) - S(t-τ))

This estimate is exact when α̂ is correct, but any error in α̂ immediately
propagates to the control output — it is a *direct* feedforward with no
smoothing.  When α varies rapidly or α̂ is noisy, Smith can produce jittery
or unstable grant rates.

IMC adds a first-order **robustness filter** on the *correction term only*:

    correction(t)   = α̂(t) · (S(t) - S(t-τ))       [unfiltered Smith term]
    corr_f(t)       = filtered version of correction(t),  τ_f = λ·T
    B̂_imc(t)        = B(t-τ) - corr_f(t)

Because the filter acts only on the correction (not on the base observation),
it does not introduce a lag in tracking B(t-τ) — it only smooths the
model-based adjustment.  When α is constant, corr_f converges to correction
and IMC gives the same result as adaptive Smith.  When α is noisy, corr_f
smooths out the noise, reducing grant-rate jitter.

**Effect of λ:**
- λ = 0: corr_f = correction instantly → same as adaptive Smith.
- λ > 0: correction is low-pass filtered with time constant λ·T.
  Larger λ = more smoothing, slower response to α changes.

State: u = [B(t), S(t), corr_f(t)]
  - B(t):      true remaining budget
  - S(t):      grant tally (as in Smith)
  - corr_f(t): low-pass filtered Smith correction term

Requires lags τ, 2τ, 3τ (for α̂ estimation).

`λ` = filter time constant as fraction of T (default 0.1).
`α` = true grant fulfilment rate (unknown to the controller).
"""
function solve_budget_imc(;
    Q     = 100.0,
    T     = 10.0,
    τ     = 1.5,
    λ     = 0.1,
    α     = 1.0,
    tspan = (0.0, T - 1e-3)
)
    τ_f = max(λ * T, τ / 20)

    function imc!(du, u, h, p, t)
        Q, T, τ, τ_f, α = p
        dt = τ > 0 ? τ : 1e-6

        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]

        # estimate α̂ from observed vs predicted drain over [t-3τ, t-2τ]
        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0

        # unfiltered Smith correction term
        correction = α_hat * (u[2] - S_now)

        # low-pass filter on correction only
        corr_f   = u[3]
        dcorr_f  = (correction - corr_f) / τ_f

        # IMC estimate: raw observation minus filtered correction
        B_hat = max(B_now - corr_f, 0.0)

        rate = B_hat / (T - t)

        du[1] = -α * rate     # true system
        du[2] =      rate     # grant tally
        du[3] = dcorr_f       # filtered correction
    end

    h(p, t) = [Q, 0.0, 0.0]   # history: full budget, zero spent, zero correction
    p = (Q, T, τ, τ_f, α)
    prob = DDEProblem(imc!, [Q, 0.0, 0.0], h, tspan, p;
                      constant_lags = [τ, 2τ, 3τ])
    solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
end
