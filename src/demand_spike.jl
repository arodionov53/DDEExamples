# ── Demand-spike variants ────────────────────────────────────────────────────
#
# Scenario: the system must spend a budget Q over horizon T, but at time
# t_spike a short burst of extra demand Δ arrives and is satisfied
# instantly (budget drops by Δ).  The spike duration d_spike < τ, so the
# delayed observation will *not* see it until t ≥ t_spike + τ.
#
# Each solver wraps the original dynamics and injects the spike as a
# one-time discontinuity in the budget state.  We model the spike by
# adjusting the initial condition for the segment [t_spike, T) while the
# history function still returns Q for t ≤ 0 and the pre-spike trajectory
# for 0 < t < t_spike.
#
# Implementation: we split the solve into two segments:
#   1. [0, t_spike) — no spike, identical to the original solver
#   2. [t_spike, T) — budget drops by Δ at t_spike; solver restarts with
#      the modified state and a history function that replays segment 1.

"""
    solve_budget_corrected_denom_spike(; Q, T, τ, t_spike, spike_Δ)

Corrected-denominator controller with a short demand spike of size `spike_Δ`
injected at `t_spike`.  The spike is instantaneous (duration → 0, which is
less than τ), so the delayed observation does not register the spike until
t ≥ t_spike + τ.
"""
function solve_budget_corrected_denom_spike(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike = 2.0, spike_Δ = 10.0
)
    # Segment 1: normal run up to t_spike
    sol1 = solve_budget_corrected_denom(; Q, T, τ, tspan = (0.0, t_spike))

    # History for segment 2: replay sol1 for t ≤ t_spike, full budget before 0
    h2(p, t) = t <= 0.0 ? [Q] : [sol1(t; idxs=1)]

    # Initial condition after spike
    B_at_spike = sol1(t_spike; idxs=1)
    u0 = [B_at_spike - spike_Δ]

    function budget_cd!(du, u, h, p, t)
        _Q, _T, _τ = p
        B_delayed = h(p, t - _τ)[1]
        du[1] = -B_delayed / (_T - t + _τ)
    end

    p = (Q, T, τ)
    prob = DDEProblem(budget_cd!, u0, h2, (t_spike, T - 1e-3), p;
                      constant_lags = [τ])
    sol2 = solve(prob, MethodOfSteps(Tsit5()))
    (sol1, sol2, t_spike)
end

"""
    solve_budget_smith_spike(; Q, T, τ, t_spike, spike_Δ)

Smith-predictor controller with a short demand spike injected at `t_spike`.
"""
function solve_budget_smith_spike(;
    Q = 100.0, T = 10.0, τ = 1.5,
    t_spike = 2.0, spike_Δ = 10.0
)
    sol1 = solve_budget_smith(; Q, T, τ, tspan = (0.0, t_spike))

    h2(p, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]

    B_at_spike = sol1(t_spike; idxs=1)
    S_at_spike = sol1(t_spike; idxs=2)
    u0 = [B_at_spike - spike_Δ, S_at_spike + spike_Δ]

    function budget_smith!(du, u, h, p, t)
        _Q, _T, _τ = p
        B_delayed = h(p, t - _τ)[1]
        S_delayed = h(p, t - _τ)[2]
        B_hat = B_delayed - (u[2] - S_delayed)
        rate  = B_hat / (_T - t)
        du[1] = -rate
        du[2] =  rate
    end

    p = (Q, T, τ)
    prob = DDEProblem(budget_smith!, u0, h2, (t_spike, T - 1e-3), p;
                      constant_lags = [τ])
    sol2 = solve(prob, MethodOfSteps(Tsit5()))
    (sol1, sol2, t_spike)
end

"""
    solve_budget_pid_pacer_spike(; Q, T, τ, t_spike, spike_Δ, Kp, Ki, Kd)

PIDPacer controller with a short demand spike injected at `t_spike`.
"""
function solve_budget_pid_pacer_spike(;
    Q = 100.0, T = 10.0, τ = 1.0,
    Kp = 1.0, Ki = 0.1, Kd = 0.05,
    max_integral = 10.0 / 0.1,
    request_rate = Q / T,
    t_spike = 2.0, spike_Δ = 10.0
)
    sol1 = solve_budget_pid_pacer(; Q, T, τ, Kp, Ki, Kd,
                                   max_integral, request_rate,
                                   tspan = (0.0, t_spike))

    h2(p, t) = t <= 0.0 ? [Q, 0.0] : [sol1(t; idxs=1), sol1(t; idxs=2)]

    target_rate = Q / T
    B_at_spike = sol1(t_spike; idxs=1)
    I_at_spike = sol1(t_spike; idxs=2)
    u0 = [B_at_spike - spike_Δ, I_at_spike]

    function pid_pacer!(du, u, h, p, t)
        _Q, _T, _τ, _Kp, _Ki, _Kd, _max_integral, _target_rate, _request_rate = p

        B_now  = h(p, t - _τ)[1]
        B_prev = h(p, t - 2_τ)[1]
        dt_obs = _τ > 0 ? _τ : 1e-6
        observed_rate = (B_prev - B_now) / dt_obs

        e = _target_rate - observed_rate
        I_new = clamp(u[2] + e, -_max_integral, _max_integral)

        B_now2 = h(p, t - 3_τ)[1]
        observed_rate_prev = (B_now - B_now2) / dt_obs
        e_prev = _target_rate - observed_rate_prev
        de = (e - e_prev) / dt_obs

        pid_out = _Kp * e + _Ki * I_new + _Kd * de
        prob    = 1.0 / (1.0 + exp(-pid_out))
        prob    = clamp(prob, 0.0, 1.0)

        du[1] = -prob * _request_rate
        du[2] = e
    end

    p = (Q, T, τ, Kp, Ki, Kd, max_integral, target_rate, request_rate)
    prob = DDEProblem(pid_pacer!, u0, h2, (t_spike, T - 1e-3), p;
                      constant_lags = [τ, 2τ, 3τ])
    sol2 = solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
    (sol1, sol2, t_spike)
end

"""
    demo_demand_spike(; Q, T, delay_fracs, t_spike_frac, spike_Δ)

Compare PIDPacer, corrected-denominator, and Smith predictor when a short
demand spike (duration < τ) hits at `t_spike_frac * T`.

Grid layout: rows = τ values (expressed as % of T), columns = controllers.
Saves `demand_spike.png`.
"""
function demo_demand_spike(;
    Q            = 100.0,
    T            = 10.0,
    delay_fracs  = [0.05, 0.10, 0.30],
    t_spike_frac = 0.2,
    spike_Δ      = 10.0,
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts      = range(0.0, T - 1e-3; length = 500)
    t_spike = t_spike_frac * T

    function stitch2(s1, s2, t_sp)
        map(t -> t < t_sp ? s1(t; idxs=1) : s2(t; idxs=1), ts)
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
        vline!(p, [t_spike];     label="spike",         linewidth=1, linestyle=:dashdot, color=:red)
        vline!(p, [t_spike + τ]; label="visible t+τ",   linewidth=1, linestyle=:dot,     color=:orange)
        plot!(p, ts, spent; linewidth=2, color=:blue,
              label="spent ($(round(spent[end]; digits=1)))")
        p
    end

    strategies = [
        ("Corr. denom",   (τ) -> begin
            s1, s2, _ = solve_budget_corrected_denom_spike(; Q, T, τ, t_spike, spike_Δ)
            stitch2(s1, s2, t_spike)
        end),
        ("Smith",         (τ) -> begin
            s1, s2, _ = solve_budget_smith_spike(; Q, T, τ, t_spike, spike_Δ)
            stitch2(s1, s2, t_spike)
        end),
        ("PID Pacer",     (τ) -> begin
            s1, s2, _ = solve_budget_pid_pacer_spike(; Q, T, τ, t_spike, spike_Δ, Kp, Ki, Kd)
            stitch2(s1, s2, t_spike)
        end),
    ]

    plts = [make_plot(name, frac * T, solver(frac * T))
            for frac in delay_fracs
            for (name, solver) in strategies]

    ncols = length(strategies)
    nrows = length(delay_fracs)
    fig = plot(plts...; layout=(nrows, ncols),
               size=(380 * ncols, 360 * nrows),
               plot_title="Demand spike (Δ=$(spike_Δ), d<τ)  Q=$(Q), T=$(T)")
    savefig(fig, "plots/demand_spike.png")
    println("Plot saved to plots/demand_spike.png")
    fig
end
