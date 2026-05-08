"""
    solve_budget_mpc_spike(; Q, T, τ, t_spike, spike_Δ)

MPC controller with a sudden budget drop of size `spike_Δ` at `t_spike`.
"""
function solve_budget_mpc_spike(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike = 2.0, spike_Δ = 10.0
)
    sol1 = solve_budget_mpc(; Q, T, τ, tspan = (0.0, t_spike))

    h2(p, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]

    B_at_spike = sol1(t_spike; idxs=1)
    S_at_spike = sol1(t_spike; idxs=2)
    u0 = [B_at_spike - spike_Δ, S_at_spike + spike_Δ]

    function mpc!(du, u, h, p, t)
        Q, T, τ = p
        dt = τ > 0 ? τ : 1e-6

        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]

        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0

        B_hat     = max(B_now - α_hat * (u[2] - S_now), 0.0)
        rate_free = B_hat  / (T - t)
        rate_safe = B_now  / (T - t)
        rate      = min(rate_free, rate_safe)

        du[1] = -rate
        du[2] =  rate
    end

    p = (Q, T, τ)
    prob = DDEProblem(mpc!, u0, h2, (t_spike, T - 1e-3), p;
                      constant_lags = [τ, 2τ, 3τ])
    sol2 = solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
    (sol1, sol2, t_spike)
end

"""
    solve_budget_imc_spike(; Q, T, τ, λ, t_spike, spike_Δ)

IMC controller with a sudden budget drop of size `spike_Δ` at `t_spike`.
"""
function solve_budget_imc_spike(;
    Q = 100.0, T = 10.0, τ = 1.5, λ = 0.1,
    t_spike = 2.0, spike_Δ = 10.0
)
    sol1 = solve_budget_imc(; Q, T, τ, λ, tspan = (0.0, t_spike))

    h2(p, t) = t <= 0.0 ? [Q, 0.0, 0.0] :
                           [sol1(t; idxs=1), sol1(t; idxs=2), sol1(t; idxs=3)]

    B_at_spike = sol1(t_spike; idxs=1)
    S_at_spike = sol1(t_spike; idxs=2)
    cf_at_spike = sol1(t_spike; idxs=3)
    u0 = [B_at_spike - spike_Δ, S_at_spike + spike_Δ, cf_at_spike]

    τ_f = max(λ * T, τ / 20)

    function imc!(du, u, h, p, t)
        Q, T, τ, τ_f = p
        dt = τ > 0 ? τ : 1e-6

        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]

        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0

        correction = α_hat * (u[2] - S_now)
        corr_f     = u[3]
        dcorr_f    = (correction - corr_f) / τ_f

        B_hat = max(B_now - corr_f, 0.0)
        rate  = B_hat / (T - t)

        du[1] = -rate
        du[2] =  rate
        du[3] = dcorr_f
    end

    p = (Q, T, τ, τ_f)
    prob = DDEProblem(imc!, u0, h2, (t_spike, T - 1e-3), p;
                      constant_lags = [τ, 2τ, 3τ])
    sol2 = solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
    (sol1, sol2, t_spike)
end

"""
    demo_mpc_spike(; Q, T, delay_fracs, t_spike_frac, spike_Δ)

Compare MPC, IMC, Smith, and corrected-denominator when a sudden budget
drop of `spike_Δ` hits at `t_spike_frac·T`.

The spike is instantaneous (duration < τ), so all controllers are blind to
it for a window of length τ.  The key question is how quickly each controller
detects and recovers from the sudden shortfall — and whether any controller
overshoots (goes below 0 or above 100% spend).

Grid layout: rows = τ values, columns = controllers.
Saves `plots/mpc_spike.png`.
"""
function demo_mpc_spike(;
    Q            = 100.0,
    T            = 10.0,
    delay_fracs  = [0.05, 0.10, 0.30],
    t_spike_frac = 0.2,
    spike_Δ      = 10.0
)
    ts      = range(0.0, T - 1e-3; length = 500)
    t_spike = t_spike_frac * T

    function stitch(s1, s2, t_sp)
        map(t -> t < t_sp ? s1(t; idxs=1) : s2(t; idxs=1), ts)
    end

    function make_plot(τ)
        τ_pct = round(Int, 100 * τ / T)

        s1_cd, s2_cd, _ = solve_budget_corrected_denom_spike(; Q, T, τ, t_spike, spike_Δ)
        s1_sm, s2_sm, _ = solve_budget_smith_spike(;          Q, T, τ, t_spike, spike_Δ)
        s1_mp, s2_mp, _ = solve_budget_mpc_spike(;            Q, T, τ, t_spike, spike_Δ)
        s1_im, s2_im, _ = solve_budget_imc_spike(;            Q, T, τ, t_spike, spike_Δ)

        B_cd = stitch(s1_cd, s2_cd, t_spike)
        B_sm = stitch(s1_sm, s2_sm, t_spike)
        B_mp = stitch(s1_mp, s2_mp, t_spike)
        B_im = stitch(s1_im, s2_im, t_spike)

        # adjust ideal line for the spike: ideal after spike accounts for Δ
        ideal = map(ts) do t
            t < t_spike ? Q * t / T :
                          Q * t / T + spike_Δ   # extra Δ was spent
        end

        p = plot(ts, ideal;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent",
            title = "τ = $(τ_pct)%·T",
            legend = :topleft)
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)
        vline!(p, [t_spike];     label = "spike −$(spike_Δ)",
               linewidth = 1, linestyle = :dashdot, color = :red)
        vline!(p, [t_spike + τ]; label = "visible t+τ",
               linewidth = 1, linestyle = :dot, color = :orange)

        for (label, B_vec, col) in [
            ("Corr. denom", B_cd, :blue),
            ("Smith",       B_sm, :green),
            ("MPC",         B_mp, :brown),
            ("IMC",         B_im, :teal),
        ]
            spent = Q .- B_vec
            plot!(p, ts, spent;
                label = "$label ($(round(spent[end]; digits=1)))",
                linewidth = 2, color = col)
        end
        p
    end

    plts = [make_plot(frac * T) for frac in delay_fracs]

    fig = plot(plts...; layout = (length(delay_fracs), 1),
               size = (800, 380 * length(delay_fracs)),
               plot_title = "Budget drop response: MPC vs IMC vs Smith vs Corr.denom  (Q=$(Q), T=$(T), Δ=$(spike_Δ))")
    savefig(fig, "plots/mpc_spike.png")
    println("Plot saved to plots/mpc_spike.png")
    fig
end
