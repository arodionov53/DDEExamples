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
    solve_budget_mpc_two_spikes(; Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)

MPC controller with two successive budget drops.
"""
function solve_budget_mpc_two_spikes(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike1 = 2.0, spike_Δ1 = 10.0,
    t_spike2 = 5.0, spike_Δ2 = 10.0
)
    function dynamics!(du, u, h, p, t)
        Q, T, τ = p
        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]
        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0
        B_hat = max(B_now - α_hat * (u[2] - S_now), 0.0)
        rate  = min(B_hat / (T - t), B_now / (T - t))
        du[1] = -rate; du[2] = rate
    end
    p = (Q, T, τ)

    sol1 = solve_budget_mpc(; Q, T, τ, tspan = (0.0, t_spike1))

    h2(p2, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]
    u0_2 = [sol1(t_spike1; idxs=1) - spike_Δ1, sol1(t_spike1; idxs=2) + spike_Δ1]
    prob2 = DDEProblem(dynamics!, u0_2, h2, (t_spike1, t_spike2), p;
                       constant_lags = [τ, 2τ, 3τ])
    sol2 = solve(prob2, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    h3(p3, t) = t <= 0.0    ? [Q, 0.0] :
                t < t_spike1 ? [sol1(t; idxs=1), sol1(t; idxs=2)] :
                               [sol2(t; idxs=1), sol2(t; idxs=2)]
    u0_3 = [sol2(t_spike2; idxs=1) - spike_Δ2, sol2(t_spike2; idxs=2) + spike_Δ2]
    prob3 = DDEProblem(dynamics!, u0_3, h3, (t_spike2, T - 1e-3), p;
                       constant_lags = [τ, 2τ, 3τ])
    sol3 = solve(prob3, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    (sol1, sol2, sol3, t_spike1, t_spike2)
end

"""
    solve_budget_imc_two_spikes(; Q, T, τ, λ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)

IMC controller with two successive budget drops.
"""
function solve_budget_imc_two_spikes(;
    Q = 100.0, T = 10.0, τ = 1.5, λ = 0.1,
    t_spike1 = 2.0, spike_Δ1 = 10.0,
    t_spike2 = 5.0, spike_Δ2 = 10.0
)
    τ_f = max(λ * T, τ / 20)

    function dynamics!(du, u, h, p, t)
        Q, T, τ, τ_f = p
        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]
        predicted_drain = S_prev - S_prev2
        observed_drain  = B_prev2 - B_prev
        α_hat   = predicted_drain > 1e-6 ?
                      clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0
        correction = α_hat * (u[2] - S_now)
        corr_f     = u[3]
        dcorr_f    = (correction - corr_f) / τ_f
        B_hat = max(B_now - corr_f, 0.0)
        rate  = B_hat / (T - t)
        du[1] = -rate; du[2] = rate; du[3] = dcorr_f
    end
    p = (Q, T, τ, τ_f)

    sol1 = solve_budget_imc(; Q, T, τ, λ, tspan = (0.0, t_spike1))

    h2(p2, t) = t <= 0.0 ? [Q, 0.0, 0.0] :
                            [sol1(t; idxs=1), sol1(t; idxs=2), sol1(t; idxs=3)]
    u0_2 = [sol1(t_spike1; idxs=1) - spike_Δ1,
            sol1(t_spike1; idxs=2) + spike_Δ1,
            sol1(t_spike1; idxs=3)]
    prob2 = DDEProblem(dynamics!, u0_2, h2, (t_spike1, t_spike2), p;
                       constant_lags = [τ, 2τ, 3τ])
    sol2 = solve(prob2, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    h3(p3, t) = t <= 0.0    ? [Q, 0.0, 0.0] :
                t < t_spike1 ? [sol1(t; idxs=1), sol1(t; idxs=2), sol1(t; idxs=3)] :
                               [sol2(t; idxs=1), sol2(t; idxs=2), sol2(t; idxs=3)]
    u0_3 = [sol2(t_spike2; idxs=1) - spike_Δ2,
            sol2(t_spike2; idxs=2) + spike_Δ2,
            sol2(t_spike2; idxs=3)]
    prob3 = DDEProblem(dynamics!, u0_3, h3, (t_spike2, T - 1e-3), p;
                       constant_lags = [τ, 2τ, 3τ])
    sol3 = solve(prob3, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    (sol1, sol2, sol3, t_spike1, t_spike2)
end

"""
    demo_mpc_spike(; Q, T, delay_fracs, t_spike_frac, spike_Δ)

Compare MPC, IMC, Smith, and corrected-denominator when a sudden budget
drop of `spike_Δ` hits at `t_spike_frac·T`.

Grid layout: one subplot per τ value.
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

    function stitch2(s1, s2, t_sp)
        map(t -> t < t_sp ? s1(t; idxs=1) : s2(t; idxs=1), ts)
    end

    function make_plot(τ)
        τ_pct = round(Int, 100 * τ / T)

        s1_cd, s2_cd, _ = solve_budget_corrected_denom_spike(; Q, T, τ, t_spike, spike_Δ)
        s1_sm, s2_sm, _ = solve_budget_smith_spike(;          Q, T, τ, t_spike, spike_Δ)
        s1_mp, s2_mp, _ = solve_budget_mpc_spike(;            Q, T, τ, t_spike, spike_Δ)
        s1_im, s2_im, _ = solve_budget_imc_spike(;            Q, T, τ, t_spike, spike_Δ)

        B_cd = stitch2(s1_cd, s2_cd, t_spike)
        B_sm = stitch2(s1_sm, s2_sm, t_spike)
        B_mp = stitch2(s1_mp, s2_mp, t_spike)
        B_im = stitch2(s1_im, s2_im, t_spike)

        ideal = map(t -> t < t_spike ? Q*t/T : Q*t/T + spike_Δ, ts)

        p = plot(ts, ideal;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent", title = "τ = $(τ_pct)%·T",
            legend = :topleft)
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)
        vline!(p, [t_spike];     label = "spike −$(spike_Δ)",
               linewidth = 1, linestyle = :dashdot, color = :red)
        vline!(p, [t_spike + τ]; label = "visible t+τ",
               linewidth = 1, linestyle = :dot, color = :orange)
        for (lbl, B_vec, col) in [("Corr. denom", B_cd, :blue),
                                   ("Smith",       B_sm, :green),
                                   ("MPC",         B_mp, :brown),
                                   ("IMC",         B_im, :teal)]
            spent = Q .- B_vec
            plot!(p, ts, spent; label = "$lbl ($(round(spent[end]; digits=1)))",
                  linewidth = 2, color = col)
        end
        p
    end

    plts = [make_plot(frac * T) for frac in delay_fracs]
    fig = plot(plts...; layout = (length(delay_fracs), 1),
               size = (800, 380 * length(delay_fracs)),
               plot_title = "Single budget drop: MPC vs IMC vs Smith vs Corr.denom  (Q=$(Q), T=$(T), Δ=$(spike_Δ))")
    savefig(fig, "plots/mpc_spike.png")
    println("Plot saved to plots/mpc_spike.png")
    fig
end

"""
    demo_mpc_two_spikes(; Q, T, delay_fracs, t_spike1_frac, spike_Δ1, t_spike2_frac, spike_Δ2)

Compare MPC, IMC, Smith, and corrected-denominator under two successive
budget drops.  Layout: one subplot per τ value.
Saves `plots/mpc_two_spikes.png`.
"""
function demo_mpc_two_spikes(;
    Q             = 100.0,
    T             = 10.0,
    delay_fracs   = [0.05, 0.10, 0.30],
    t_spike1_frac = 0.2,
    spike_Δ1      = 10.0,
    t_spike2_frac = 0.5,
    spike_Δ2      = 10.0
)
    ts       = range(0.0, T - 1e-3; length = 500)
    t_spike1 = t_spike1_frac * T
    t_spike2 = t_spike2_frac * T

    function stitch3(s1, s2, s3, ts1, ts2)
        map(t -> t < ts1 ? s1(t; idxs=1) :
                 t < ts2 ? s2(t; idxs=1) :
                           s3(t; idxs=1), ts)
    end

    function make_plot(τ)
        τ_pct = round(Int, 100 * τ / T)

        s1_cd, s2_cd, s3_cd, _, _ = solve_budget_corrected_denom_two_spikes(;
            Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)
        s1_sm, s2_sm, s3_sm, _, _ = solve_budget_smith_two_spikes(;
            Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)
        s1_mp, s2_mp, s3_mp, _, _ = solve_budget_mpc_two_spikes(;
            Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)
        s1_im, s2_im, s3_im, _, _ = solve_budget_imc_two_spikes(;
            Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)

        B_cd = stitch3(s1_cd, s2_cd, s3_cd, t_spike1, t_spike2)
        B_sm = stitch3(s1_sm, s2_sm, s3_sm, t_spike1, t_spike2)
        B_mp = stitch3(s1_mp, s2_mp, s3_mp, t_spike1, t_spike2)
        B_im = stitch3(s1_im, s2_im, s3_im, t_spike1, t_spike2)

        total_Δ = spike_Δ1 + spike_Δ2
        ideal = map(ts) do t
            t < t_spike1 ? Q*t/T :
            t < t_spike2 ? Q*t/T + spike_Δ1 :
                           Q*t/T + total_Δ
        end

        p = plot(ts, ideal;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent", title = "τ = $(τ_pct)%·T",
            legend = :topleft)
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)
        vline!(p, [t_spike1];       label = "spike1 −$(spike_Δ1)",
               linewidth = 1, linestyle = :dashdot, color = :red)
        vline!(p, [t_spike1 + τ];   label = "spike1 vis. t+τ",
               linewidth = 1, linestyle = :dot,     color = :orange)
        vline!(p, [t_spike2];       label = "spike2 −$(spike_Δ2)",
               linewidth = 1, linestyle = :dashdot, color = :darkred)
        vline!(p, [t_spike2 + τ];   label = "spike2 vis. t+τ",
               linewidth = 1, linestyle = :dot,     color = :darkorange)
        for (lbl, B_vec, col) in [("Corr. denom", B_cd, :blue),
                                   ("Smith",       B_sm, :green),
                                   ("MPC",         B_mp, :brown),
                                   ("IMC",         B_im, :teal)]
            spent = Q .- B_vec
            plot!(p, ts, spent; label = "$lbl ($(round(spent[end]; digits=1)))",
                  linewidth = 2, color = col)
        end
        p
    end

    plts = [make_plot(frac * T) for frac in delay_fracs]
    fig = plot(plts...; layout = (length(delay_fracs), 1),
               size = (800, 400 * length(delay_fracs)),
               plot_title = "Two budget drops: MPC vs IMC vs Smith vs Corr.denom  (Q=$(Q), T=$(T), Δ1=$(spike_Δ1), Δ2=$(spike_Δ2))")
    savefig(fig, "plots/mpc_two_spikes.png")
    println("Plot saved to plots/mpc_two_spikes.png")
    fig
end
