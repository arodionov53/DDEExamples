# ── Two-spike variants ───────────────────────────────────────────────────────
#
# Extension of the single-spike pattern to two successive instantaneous
# demand spikes, both with duration < τ so neither is visible in the delayed
# observation before it has already passed.
#
# Implementation: three segments separated at t_spike1 < t_spike2:
#   seg1  [0,        t_spike1)  — normal dynamics
#   seg2  [t_spike1, t_spike2)  — budget drops by Δ1; history replays seg1
#   seg3  [t_spike2, T)         — budget drops by Δ2; history replays seg1+seg2

"""
    solve_budget_corrected_denom_two_spikes(; Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)

Corrected-denominator controller with two short demand spikes.
"""
function solve_budget_corrected_denom_two_spikes(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike1 = 2.0, spike_Δ1 = 10.0,
    t_spike2 = 5.0, spike_Δ2 = 10.0
)
    function dynamics!(du, u, h, p, t)
        _Q, _T, _τ = p
        du[1] = -h(p, t - _τ)[1] / (_T - t + _τ)
    end
    p = (Q, T, τ)

    sol1 = solve_budget_corrected_denom(; Q, T, τ, tspan = (0.0, t_spike1))

    h2(p2, t) = t <= 0.0 ? [Q] : [sol1(t; idxs=1)]
    u0_2 = [sol1(t_spike1; idxs=1) - spike_Δ1]
    prob2 = DDEProblem(dynamics!, u0_2, h2, (t_spike1, t_spike2), p; constant_lags=[τ])
    sol2 = solve(prob2, MethodOfSteps(Tsit5()))

    h3(p3, t) = t <= 0.0   ? [Q] :
                t < t_spike1 ? [sol1(t; idxs=1)] :
                               [sol2(t; idxs=1)]
    u0_3 = [sol2(t_spike2; idxs=1) - spike_Δ2]
    prob3 = DDEProblem(dynamics!, u0_3, h3, (t_spike2, T - 1e-3), p; constant_lags=[τ])
    sol3 = solve(prob3, MethodOfSteps(Tsit5()))

    (sol1, sol2, sol3, t_spike1, t_spike2)
end

"""
    solve_budget_smith_two_spikes(; Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)

Smith-predictor controller with two short demand spikes.
"""
function solve_budget_smith_two_spikes(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike1 = 2.0, spike_Δ1 = 10.0,
    t_spike2 = 5.0, spike_Δ2 = 10.0
)
    function dynamics!(du, u, h, p, t)
        _Q, _T, _τ = p
        B_del = h(p, t - _τ)[1]
        S_del = h(p, t - _τ)[2]
        B_hat = B_del - (u[2] - S_del)
        rate  = B_hat / (_T - t)
        du[1] = -rate; du[2] = rate
    end
    p = (Q, T, τ)

    sol1 = solve_budget_smith(; Q, T, τ, tspan = (0.0, t_spike1))

    h2(p2, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]
    u0_2 = [sol1(t_spike1; idxs=1) - spike_Δ1, sol1(t_spike1; idxs=2) + spike_Δ1]
    prob2 = DDEProblem(dynamics!, u0_2, h2, (t_spike1, t_spike2), p; constant_lags=[τ])
    sol2 = solve(prob2, MethodOfSteps(Tsit5()))

    h3(p3, t) = t <= 0.0   ? [Q, 0.0] :
                t < t_spike1 ? [sol1(t; idxs=1), sol1(t; idxs=2)] :
                               [sol2(t; idxs=1), sol2(t; idxs=2)]
    u0_3 = [sol2(t_spike2; idxs=1) - spike_Δ2, sol2(t_spike2; idxs=2) + spike_Δ2]
    prob3 = DDEProblem(dynamics!, u0_3, h3, (t_spike2, T - 1e-3), p; constant_lags=[τ])
    sol3 = solve(prob3, MethodOfSteps(Tsit5()))

    (sol1, sol2, sol3, t_spike1, t_spike2)
end

"""
    solve_budget_pid_pacer_two_spikes(; Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2, Kp, Ki, Kd)

PIDPacer controller with two short demand spikes.
"""
function solve_budget_pid_pacer_two_spikes(;
    Q = 100.0, T = 10.0, τ = 1.0,
    Kp = 1.0, Ki = 0.1, Kd = 0.05,
    max_integral = 10.0 / 0.1,
    request_rate = Q / T,
    t_spike1 = 2.0, spike_Δ1 = 10.0,
    t_spike2 = 5.0, spike_Δ2 = 10.0
)
    target_rate = Q / T

    function dynamics!(du, u, h, p, t)
        _Q, _T, _τ, _Kp, _Ki, _Kd, _max_integral, _target_rate, _request_rate = p
        dt_obs = _τ > 0 ? _τ : 1e-6
        B_now  = h(p, t -   _τ)[1]
        B_prev = h(p, t - 2_τ)[1]
        B_prev2 = h(p, t - 3_τ)[1]
        observed_rate      = (B_prev  - B_now ) / dt_obs
        observed_rate_prev = (B_prev2 - B_prev) / dt_obs
        e      = _target_rate - observed_rate
        e_prev = _target_rate - observed_rate_prev
        I_new  = clamp(u[2] + e, -_max_integral, _max_integral)
        de     = (e - e_prev) / dt_obs
        pid_out = _Kp * e + _Ki * I_new + _Kd * de
        prob    = clamp(1.0 / (1.0 + exp(-pid_out)), 0.0, 1.0)
        du[1]   = -prob * _request_rate
        du[2]   = e
    end
    p = (Q, T, τ, Kp, Ki, Kd, max_integral, target_rate, request_rate)

    sol1 = solve_budget_pid_pacer(; Q, T, τ, Kp, Ki, Kd,
                                   max_integral, request_rate,
                                   tspan = (0.0, t_spike1))

    h2(p2, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]
    u0_2 = [sol1(t_spike1; idxs=1) - spike_Δ1, sol1(t_spike1; idxs=2)]
    prob2 = DDEProblem(dynamics!, u0_2, h2, (t_spike1, t_spike2), p;
                       constant_lags=[τ, 2τ, 3τ])
    sol2 = solve(prob2, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    h3(p3, t) = t <= 0.0   ? [Q, 0.0] :
                t < t_spike1 ? [sol1(t; idxs=1), sol1(t; idxs=2)] :
                               [sol2(t; idxs=1), sol2(t; idxs=2)]
    u0_3 = [sol2(t_spike2; idxs=1) - spike_Δ2, sol2(t_spike2; idxs=2)]
    prob3 = DDEProblem(dynamics!, u0_3, h3, (t_spike2, T - 1e-3), p;
                       constant_lags=[τ, 2τ, 3τ])
    sol3 = solve(prob3, MethodOfSteps(Tsit5()); dtmax = τ / 10)

    (sol1, sol2, sol3, t_spike1, t_spike2)
end

"""
    demo_demand_two_spikes(; Q, T, delay_fracs, t_spike1_frac, spike_Δ1, t_spike2_frac, spike_Δ2)

Compare PIDPacer, corrected-denominator, and Smith predictor under two
successive instantaneous demand spikes, each shorter than τ.

The second spike may arrive before the controller has fully recovered from
the first (when t_spike2 - t_spike1 < τ), stressing the blind-window
response further.

Grid layout: rows = τ values (expressed as % of T), columns = controllers.
Saves `demand_two_spikes.png`.
"""
function demo_demand_two_spikes(;
    Q              = 100.0,
    T              = 10.0,
    delay_fracs    = [0.05, 0.10, 0.30],
    t_spike1_frac  = 0.2,
    spike_Δ1       = 10.0,
    t_spike2_frac  = 0.5,
    spike_Δ2       = 10.0,
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts      = range(0.0, T - 1e-3; length = 500)
    t_spike1 = t_spike1_frac * T
    t_spike2 = t_spike2_frac * T

    function stitch3(s1, s2, s3, ts1, ts2)
        map(t -> t < ts1 ? s1(t; idxs=1) :
                 t < ts2 ? s2(t; idxs=1) :
                           s3(t; idxs=1), ts)
    end

    function make_plot(controller, τ, B_vec)
        τ_pct = round(Int, 100 * τ / T)
        spent = Q .- B_vec
        p = plot(ts, Q .* ts ./ T;
            label="ideal", linewidth=2, linestyle=:dash, color=:black,
            xlabel="time", ylabel="spent",
            title="$(controller)  τ=$(τ_pct)%·T",
            legend=:topleft)
        hline!(p, [Q]; label="cap", linewidth=1, linestyle=:dot, color=:grey)
        vline!(p, [t_spike1];       label="spike1",           linewidth=1, linestyle=:dashdot, color=:red)
        vline!(p, [t_spike1 + τ];   label="spike1 vis. t+τ",  linewidth=1, linestyle=:dot,     color=:orange)
        vline!(p, [t_spike2];       label="spike2",           linewidth=1, linestyle=:dashdot, color=:darkred)
        vline!(p, [t_spike2 + τ];   label="spike2 vis. t+τ",  linewidth=1, linestyle=:dot,     color=:darkorange)
        plot!(p, ts, spent; linewidth=2, color=:blue,
              label="spent ($(round(spent[end]; digits=1)))")
        p
    end

    strategies = [
        ("Corr. denom", (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_corrected_denom_two_spikes(;
                Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
        ("Smith",       (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_smith_two_spikes(;
                Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
        ("PID Pacer",   (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_pid_pacer_two_spikes(;
                Q, T, τ, t_spike1, spike_Δ1, t_spike2, spike_Δ2, Kp, Ki, Kd)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
    ]

    plts = [make_plot(name, frac * T, solver(frac * T))
            for frac in delay_fracs
            for (name, solver) in strategies]

    ncols = length(strategies)
    nrows = length(delay_fracs)
    fig = plot(plts...; layout=(nrows, ncols),
               size=(380 * ncols, 360 * nrows),
               plot_title="Two demand spikes (Δ1=$(spike_Δ1), Δ2=$(spike_Δ2), d<τ)  Q=$(Q), T=$(T)")
    savefig(fig, "demand_two_spikes.png")
    println("Plot saved to demand_two_spikes.png")
    fig
end

"""
    demo_demand_spike_then_drop(; Q, T, delay_fracs, t_spike_frac, spike_Δ, t_drop_frac, drop_Δ)

Like `demo_demand_two_spikes` but the second event is a sudden budget *refund*
(drop in demand) of size `drop_Δ` — modelled as a negative spike, i.e. `B`
increases by `drop_Δ` instantaneously at `t_drop`.

This tests whether controllers can exploit a windfall: level-based controllers
(Smith, corrected denom) should slow spending immediately; the rate-based
PIDPacer will mis-read the refund as a sudden deceleration and *increase*
grant probability, potentially over-spending.

Grid layout: rows = τ values (expressed as % of T), columns = controllers.
Saves `demand_spike_then_drop.png`.
"""
function demo_demand_spike_then_drop(;
    Q             = 100.0,
    T             = 10.0,
    delay_fracs   = [0.05, 0.10, 0.30],
    t_spike_frac  = 0.2,
    spike_Δ       = 10.0,
    t_drop_frac   = 0.5,
    drop_Δ        = 10.0,    # positive value → budget refunded
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts      = range(0.0, T - 1e-3; length = 500)
    t_spike = t_spike_frac * T
    t_drop  = t_drop_frac  * T

    function stitch3(s1, s2, s3, ts1, ts2)
        map(t -> t < ts1 ? s1(t; idxs=1) :
                 t < ts2 ? s2(t; idxs=1) :
                           s3(t; idxs=1), ts)
    end

    function make_plot(controller, τ, B_vec)
        τ_pct = round(Int, 100 * τ / T)
        spent = Q .- B_vec
        p = plot(ts, Q .* ts ./ T;
            label="ideal", linewidth=2, linestyle=:dash, color=:black,
            xlabel="time", ylabel="spent",
            title="$(controller)  τ=$(τ_pct)%·T",
            legend=:topleft)
        hline!(p, [Q]; label="cap", linewidth=1, linestyle=:dot, color=:grey)
        vline!(p, [t_spike];        label="spike (−$(spike_Δ))",  linewidth=1, linestyle=:dashdot, color=:red)
        vline!(p, [t_spike + τ];    label="spike vis. t+τ",       linewidth=1, linestyle=:dot,     color=:orange)
        vline!(p, [t_drop];         label="drop (+$(drop_Δ))",    linewidth=1, linestyle=:dashdot, color=:darkgreen)
        vline!(p, [t_drop + τ];     label="drop vis. t+τ",        linewidth=1, linestyle=:dot,     color=:green)
        plot!(p, ts, spent; linewidth=2, color=:blue,
              label="spent ($(round(spent[end]; digits=1)))")
        p
    end

    # Reuse two-spike solvers: spike_Δ1 = +spike (budget down),
    # spike_Δ2 = -drop_Δ (budget up, refund).
    strategies = [
        ("Corr. denom", (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_corrected_denom_two_spikes(;
                Q, T, τ, t_spike1=t_spike, spike_Δ1=spike_Δ,
                         t_spike2=t_drop,  spike_Δ2=-drop_Δ)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
        ("Smith",       (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_smith_two_spikes(;
                Q, T, τ, t_spike1=t_spike, spike_Δ1=spike_Δ,
                         t_spike2=t_drop,  spike_Δ2=-drop_Δ)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
        ("PID Pacer",   (τ) -> begin
            s1,s2,s3,ts1,ts2 = solve_budget_pid_pacer_two_spikes(;
                Q, T, τ, t_spike1=t_spike, spike_Δ1=spike_Δ,
                         t_spike2=t_drop,  spike_Δ2=-drop_Δ, Kp, Ki, Kd)
            stitch3(s1, s2, s3, ts1, ts2)
        end),
    ]

    plts = [make_plot(name, frac * T, solver(frac * T))
            for frac in delay_fracs
            for (name, solver) in strategies]

    ncols = length(strategies)
    nrows = length(delay_fracs)
    fig = plot(plts...; layout=(nrows, ncols),
               size=(380 * ncols, 360 * nrows),
               plot_title="Spike then refund (spike −$(spike_Δ) @ $(t_spike_frac)T, drop +$(drop_Δ) @ $(t_drop_frac)T)  Q=$(Q), T=$(T)")
    savefig(fig, "demand_spike_then_drop.png")
    println("Plot saved to demand_spike_then_drop.png")
    fig
end
