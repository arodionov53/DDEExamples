# ── Random-demand helpers ─────────────────────────────────────────────────────
#
# Solve a budget-pacing DDE under a sequence of N instantaneous demand events
# (spikes and drops) at times t_events with signed sizes Δ_events (positive =
# budget withdrawn, negative = budget refunded).
#
# Each event is a discontinuity; we solve segment-by-segment between them.
# Returns a vector of segment solutions and the event times.

function _solve_random_demand_cd(; Q, T, τ, t_events, Δ_events)
    function dyn!(du, u, h, p, t)
        _Q, _T, _τ = p
        du[1] = -h(p, t - _τ)[1] / (_T - t + _τ)
    end
    p = (Q, T, τ)
    segs = []
    boundaries = [0.0; t_events; T - 1e-3]

    sol_prev = solve_budget_corrected_denom(; Q, T, τ, tspan=(0.0, boundaries[2]))
    push!(segs, sol_prev)

    for k in 2:length(boundaries)-1
        t0, t1 = boundaries[k], boundaries[k+1]
        idx = k - 1   # index into t_events / Δ_events
        segs_so_far = copy(segs)
        # history: replay all prior segments
        function h_k(pp, t)
            t <= 0.0 && return [Q]
            for (i, seg) in enumerate(segs_so_far)
                if t < boundaries[i+1]
                    return [seg(t; idxs=1)]
                end
            end
            return [segs_so_far[end](t; idxs=1)]
        end
        B0 = segs[end](t0; idxs=1) - Δ_events[idx]
        prob = DDEProblem(dyn!, [B0], h_k, (t0, t1), p; constant_lags=[τ])
        sol = solve(prob, MethodOfSteps(Tsit5()))
        sol.retcode == ReturnCode.Success || return (segs, boundaries)
        push!(segs, sol)
    end
    (segs, boundaries)
end

function _solve_random_demand_smith(; Q, T, τ, t_events, Δ_events)
    function dyn!(du, u, h, p, t)
        _Q, _T, _τ = p
        B_del = h(p, t - _τ)[1]
        S_del = h(p, t - _τ)[2]
        B_hat = B_del - (u[2] - S_del)
        rate  = B_hat / (_T - t)
        du[1] = -rate; du[2] = rate
    end
    p = (Q, T, τ)
    segs = []
    boundaries = [0.0; t_events; T - 1e-3]

    sol_prev = solve_budget_smith(; Q, T, τ, tspan=(0.0, boundaries[2]))
    push!(segs, sol_prev)

    for k in 2:length(boundaries)-1
        t0, t1 = boundaries[k], boundaries[k+1]
        idx = k - 1
        segs_so_far = copy(segs)
        function h_k(pp, t)
            t <= 0.0 && return [Q, 0.0]
            for (i, seg) in enumerate(segs_so_far)
                if t < boundaries[i+1]
                    return [seg(t; idxs=1), seg(t; idxs=2)]
                end
            end
            return [segs_so_far[end](t; idxs=1), segs_so_far[end](t; idxs=2)]
        end
        B0 = segs[end](t0; idxs=1) - Δ_events[idx]
        S0 = segs[end](t0; idxs=2) + Δ_events[idx]
        prob = DDEProblem(dyn!, [B0, S0], h_k, (t0, t1), p; constant_lags=[τ])
        sol = solve(prob, MethodOfSteps(Tsit5()))
        sol.retcode == ReturnCode.Success || return (segs, boundaries)
        push!(segs, sol)
    end
    (segs, boundaries)
end

function _solve_random_demand_pid_pacer(; Q, T, τ, t_events, Δ_events,
    Kp=1.0, Ki=0.1, Kd=0.05, max_integral=10.0/0.1, request_rate=Q/T)
    target_rate = Q / T
    function dyn!(du, u, h, p, t)
        _Q, _T, _τ, _Kp, _Ki, _Kd, _max_integral, _target_rate, _request_rate = p
        dt_obs = _τ > 0 ? _τ : 1e-6
        B_now   = h(p, t -   _τ)[1]
        B_prev  = h(p, t - 2_τ)[1]
        B_prev2 = h(p, t - 3_τ)[1]
        obs_rate      = (B_prev  - B_now ) / dt_obs
        obs_rate_prev = (B_prev2 - B_prev) / dt_obs
        e      = _target_rate - obs_rate
        e_prev = _target_rate - obs_rate_prev
        I_new  = clamp(u[2] + e, -_max_integral, _max_integral)
        de     = (e - e_prev) / dt_obs
        pid_out = _Kp * e + _Ki * I_new + _Kd * de
        prob    = clamp(1.0 / (1.0 + exp(-pid_out)), 0.0, 1.0)
        du[1] = -prob * _request_rate
        du[2] = e
    end
    p = (Q, T, τ, Kp, Ki, Kd, max_integral, target_rate, request_rate)
    segs = []
    boundaries = [0.0; t_events; T - 1e-3]

    sol_prev = solve_budget_pid_pacer(; Q, T, τ, Kp, Ki, Kd,
                                       max_integral, request_rate,
                                       tspan=(0.0, boundaries[2]))
    push!(segs, sol_prev)

    for k in 2:length(boundaries)-1
        t0, t1 = boundaries[k], boundaries[k+1]
        idx = k - 1
        segs_so_far = copy(segs)
        function h_k(pp, t)
            t <= 0.0 && return [Q, 0.0]
            for (i, seg) in enumerate(segs_so_far)
                if t < boundaries[i+1]
                    return [seg(t; idxs=1), seg(t; idxs=2)]
                end
            end
            return [segs_so_far[end](t; idxs=1), segs_so_far[end](t; idxs=2)]
        end
        B0 = segs[end](t0; idxs=1) - Δ_events[idx]
        I0 = segs[end](t0; idxs=2)
        prob = DDEProblem(dyn!, [B0, I0], h_k, (t0, t1), p;
                          constant_lags=[τ, 2τ, 3τ])
        sol = solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
        sol.retcode == ReturnCode.Success || return (segs, boundaries)
        push!(segs, sol)
    end
    (segs, boundaries)
end

# Stitch N segment solutions into a single vector over evaluation times ts.
function _stitch_segs(segs, boundaries, ts; idxs=1)
    map(ts) do t
        for (i, seg) in enumerate(segs)
            if i == length(segs) || t < boundaries[i+1]
                return seg(t; idxs=idxs)
            end
        end
        segs[end](t; idxs=idxs)
    end
end

"""
    demo_demand_random(; Q, T, delay_fracs, n_events, event_Δ_max, n_samples)

Compare PIDPacer, corrected-denominator, and Smith predictor under a random
sequence of instantaneous demand events (spikes and drops) whose arrival times
and signs are drawn independently for each Monte Carlo trial.

Each trial draws:
  - `n_events` arrival times ~ Uniform(0, T), sorted.
  - Each event size ~ Uniform(0, event_Δ_max) with random sign (±),
    clamped so B never goes negative.

The solve is split at every event time (N+1 segments).  `n_samples` trials are
run per controller per τ; the plot shows the mean trajectory ± 1σ ribbon.

Grid layout: rows = τ values (% of T), columns = controllers.
Saves `demand_random.png`.
"""
function demo_demand_random(;
    Q            = 100.0,
    T            = 10.0,
    delay_fracs  = [0.05, 0.10, 0.30],
    n_events     = 5,
    event_Δ_max  = 8.0,
    n_samples    = 30,
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts = range(0.0, T - 1e-3; length = 500)

    strategies = [
        ("Corr. denom", (τ, t_ev, Δ_ev) -> begin
            segs, bounds = _solve_random_demand_cd(;
                Q, T, τ, t_events=t_ev, Δ_events=Δ_ev)
            _stitch_segs(segs, bounds, ts)
        end),
        ("Smith",       (τ, t_ev, Δ_ev) -> begin
            segs, bounds = _solve_random_demand_smith(;
                Q, T, τ, t_events=t_ev, Δ_events=Δ_ev)
            _stitch_segs(segs, bounds, ts)
        end),
        ("PID Pacer",   (τ, t_ev, Δ_ev) -> begin
            segs, bounds = _solve_random_demand_pid_pacer(;
                Q, T, τ, t_events=t_ev, Δ_events=Δ_ev, Kp, Ki, Kd)
            _stitch_segs(segs, bounds, ts)
        end),
    ]

    plts = []
    for frac in delay_fracs
        τ = frac * T
        τ_pct = round(Int, 100 * frac)
        for (name, solver) in strategies
            cols = Vector{Vector{Float64}}()
            while length(cols) < n_samples
                # Draw random event times in [0, 0.7·T]; keeping events
                # away from the deadline avoids B̂/(T-t) stiffening.
                t_max_event = T * 0.70
                t_ev = sort(rand(n_events) .* t_max_event)
                # Draw signed sizes; total withdrawals capped at 40% of Q
                # so the solver never faces a near-zero budget near the deadline.
                signs = 2 .* (rand(n_events) .> 0.5) .- 1
                sizes = rand(n_events) .* event_Δ_max
                Δ_ev  = signs .* sizes
                sum(Δ_ev) > 0.4 * Q && continue   # net withdrawal too large, resample
                try
                    B_vec = solver(τ, t_ev, Δ_ev)
                    sol_ok = all(isfinite, B_vec) && minimum(B_vec) > 0.0
                    sol_ok && push!(cols, Q .- B_vec)
                catch
                end
            end
            mat = hcat(cols...)
            μ = vec(mean(mat; dims=2))
            σ = vec(std(mat;  dims=2))
            spent_end = round(μ[end]; digits=1)

            p = plot(ts, Q .* ts ./ T;
                label="ideal", linewidth=2, linestyle=:dash, color=:black,
                xlabel="time", ylabel="spent",
                title="$(name)  τ=$(τ_pct)%·T",
                legend=:topleft)
            hline!(p, [Q]; label="cap", linewidth=1, linestyle=:dot, color=:grey)
            plot!(p, ts, μ; ribbon=σ, fillalpha=0.2, linewidth=2, color=:blue,
                  label="mean±σ ($(spent_end))")
            push!(plts, p)
        end
    end

    ncols = length(strategies)
    nrows = length(delay_fracs)
    fig = plot(plts...; layout=(nrows, ncols),
               size=(380 * ncols, 360 * nrows),
               plot_title="Random demand events (n=$(n_events), Δ≤$(event_Δ_max))  Q=$(Q), T=$(T)")
    savefig(fig, "plots/demand_random.png")
    println("Plot saved to plots/demand_random.png")
    fig
end
