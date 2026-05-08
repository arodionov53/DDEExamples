"""
    solve_budget_delay(; Q, T, τ, tspan)

Budget spending with information delay.

A controller must spend a total budget Q evenly over a horizon T, but only
observes the remaining budget with a delay τ.  At each moment the controller
targets the rate that would exhaust the *observed* remaining budget by time T:

    dB/dt = -B(t-τ) / (T - t)

History: B(t) = Q for t ≤ 0 (full budget available at start).

* τ = 0  → perfect pacing, B(t) = Q·(1 - t/T).
* τ > 0  → the controller always "thinks" it has more than it does, so it
           under-spends early and then rushes near the deadline.

Returns the DDE solution.  Call `demo_budget_delay()` to see the comparison plot.
"""
function solve_budget_delay(;
    Q = 100.0, T = 10.0, τ = 1.5,
    tspan = (0.0, T - 1e-3)   # stop just before T to avoid division by zero
)
    function budget!(du, u, h, p, t)
        Q, T, τ = p
        B_delayed = h(p, t - τ)[1]
        remaining_time = T - t
        du[1] = -B_delayed / remaining_time
    end

    h(p, t) = [Q]          # full budget in the past
    p = (Q, T, τ)
    prob = DDEProblem(budget!, [Q], h, tspan, p; constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

"""
    solve_budget_nodelay(; Q, T)

Ideal budget spending with zero delay (ODE): dB/dt = -B(t) / (T - t).

When τ = 0 the delayed balance B(t-τ) equals B(t) and the DDE reduces to
this separable ODE with exact solution B(t) = Q·(1 - t/T) — perfect linear
drawdown that reaches zero exactly at the deadline.
"""
function solve_budget_nodelay(;
    Q = 100.0, T = 10.0,
    tspan = (0.0, T - 1e-3)
)
    function budget_ode!(du, u, p, t)
        Q, T = p
        du[1] = -u[1] / (T - t)
    end

    p = (Q, T)
    prob = ODEProblem(budget_ode!, [Q], tspan, p)
    solve(prob, Tsit5())
end

"""
    demo_budget_nodelay(; Q, T, delays)

Plot the τ = 0 ideal ODE solution alongside several naive-DDE trajectories.

Two subplots:
  - Left:  remaining balance B(t) — ideal straight line vs. delayed curves
  - Right: cumulative spent Q - B(t) — shows how delay causes late-stage rush

The exact analytical solution B*(t) = Q·(1 - t/T) is drawn as a thin dotted
line on top of the ODE solution to confirm they coincide.

Saves `plots/budget_nodelay.png`.
"""
function demo_budget_nodelay(;
    Q = 100.0, T = 10.0,
    delays = [0.05, 0.10, 0.30] .* T
)
    ts = range(0.0, T - 1e-3; length = 500)
    sol0 = solve_budget_nodelay(; Q, T)
    B0   = sol0.(ts; idxs = 1)
    exact = Q .* (1 .- ts ./ T)

    colors = [:red, :orange, :purple]

    # ── left subplot: balance ─────────────────────────────────────────────
    p1 = plot(ts, B0;
        label = "τ = 0 (ODE)", linewidth = 2, color = :blue,
        xlabel = "time", ylabel = "remaining budget B(t)",
        title = "Remaining budget",
        legend = :topright)
    plot!(p1, ts, exact;
        label = "exact B*(t) = Q(1-t/T)", linewidth = 1,
        linestyle = :dot, color = :black)
    hline!(p1, [0.0]; label = "", linewidth = 1, linestyle = :dash, color = :grey)

    for (τ, col) in zip(delays, colors)
        sol = solve_budget_delay(; Q, T, τ)
        τ_pct = round(Int, 100 * τ / T)
        B = sol.(ts; idxs = 1)
        plot!(p1, ts, B; label = "naive τ=$(τ_pct)%·T", linewidth = 2, color = col)
    end

    # ── right subplot: spent ──────────────────────────────────────────────
    p2 = plot(ts, Q .- B0;
        label = "τ = 0 (ODE)", linewidth = 2, color = :blue,
        xlabel = "time", ylabel = "cumulative spent Q - B(t)",
        title = "Cumulative spend",
        legend = :topleft)
    plot!(p2, ts, Q .* ts ./ T;
        label = "ideal linear", linewidth = 1, linestyle = :dot, color = :black)
    hline!(p2, [Q]; label = "budget cap", linewidth = 1, linestyle = :dash, color = :grey)

    for (τ, col) in zip(delays, colors)
        sol = solve_budget_delay(; Q, T, τ)
        τ_pct = round(Int, 100 * τ / T)
        spent = Q .- sol.(ts; idxs = 1)
        plot!(p2, ts, spent; label = "naive τ=$(τ_pct)%·T", linewidth = 2, color = col)
    end

    fig = plot(p1, p2; layout = (1, 2), size = (1000, 420),
               plot_title = "τ = 0 ideal vs. naive delayed controller  (Q=$Q, T=$T)")
    savefig(fig, "plots/budget_nodelay.png")
    println("Plot saved to plots/budget_nodelay.png")
    fig
end

"""
    solve_budget_corrected_denom(; Q, T, τ)

Corrected-denominator controller: divide by (T - t + τ) instead of (T - t).

    dB/dt = -B(t-τ) / (T - t + τ)

The extra τ in the denominator compensates for the fact that the observed
balance is τ time units stale.  As τ → 0 it reduces to the ideal controller.
"""
function solve_budget_corrected_denom(;
    Q = 100.0, T = 10.0, τ = 1.5,
    tspan = (0.0, T - 1e-3)
)
    function budget_cd!(du, u, h, p, t)
        Q, T, τ = p
        B_delayed = h(p, t - τ)[1]
        du[1] = -B_delayed / (T - t + τ)
    end

    h(p, t) = [Q]
    p = (Q, T, τ)
    prob = DDEProblem(budget_cd!, [Q], h, tspan, p; constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

"""
    solve_budget_smith(; Q, T, τ)

Smith-predictor controller: reconstruct the true current budget from the
delayed observation and the known cumulative spend since t-τ.

State: u = [B(t), S(t)] where S(t) = total spent up to t = Q - B(t).

    B̂(t) = B(t-τ) - (S(t) - S(t-τ))   ← predicted current balance
    dB/dt = -B̂(t) / (T - t)
    dS/dt = -dB/dt

Since B̂(t) = B(t) exactly (the prediction is perfect when the model is
correct), this reduces to the ideal ODE and B(T) = 0.
"""
function solve_budget_smith(;
    Q = 100.0, T = 10.0, τ = 1.5,
    tspan = (0.0, T - 1e-3)
)
    function budget_smith!(du, u, h, p, t)
        Q, T, τ = p
        B_delayed = h(p, t - τ)[1]
        S_delayed = h(p, t - τ)[2]
        B_hat = B_delayed - (u[2] - S_delayed)   # predicted current balance
        rate = B_hat / (T - t)
        du[1] = -rate
        du[2] =  rate
    end

    h(p, t) = [Q, 0.0]   # full budget, zero spent in history
    p = (Q, T, τ)
    prob = DDEProblem(budget_smith!, [Q, 0.0], h, tspan, p; constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

"""
    solve_budget_pid(; Q, T, τ, Kp, Ki, Kd)

PID controller for budget spending under information delay.

State: u = [B(t), I(t)] where I(t) is the integral of the tracking error.

The reference trajectory is B_ref(t) = Q·(1 - t/T) (ideal linear drawdown).
The observed error uses the delayed balance:

    e(t)     = B(t-τ) - B_ref(t-τ)        # delayed tracking error
    spend(t) = Q/T + Kp·e(t) + Ki·I(t) + Kd·ė(t)

where ė(t) = (e(t) - e(t-τ)) / τ is approximated via a second delay.
I(t) integrates the delayed error to remove steady-state offset.

Returns the DDE solution (index 1 = B, index 2 = I).
"""
function solve_budget_pid(;
    Q = 100.0, T = 10.0, τ = 1.5,
    Kp = 1.0, Ki = 0.5, Kd = 0.1,
    tspan = (0.0, T - 1e-3)
)
    function pid!(du, u, h, p, t)
        Q, T, τ, Kp, Ki, Kd = p
        B_ref(s) = Q * (1 - s / T)

        B_del  = h(p, t - τ)[1]
        B_del2 = h(p, t - 2τ)[1]   # one extra delay for derivative estimate

        e      = B_del  - B_ref(t - τ)
        e_prev = B_del2 - B_ref(t - 2τ)
        de     = τ > 0 ? (e - e_prev) / τ : 0.0

        spend = Q / T + Kp * e + Ki * u[2] + Kd * de
        du[1] = -spend
        du[2] = e
    end

    h(p, t) = [Q, 0.0]
    p = (Q, T, τ, Kp, Ki, Kd)
    prob = DDEProblem(pid!, [Q, 0.0], h, tspan, p; constant_lags = [τ, 2.0*τ])
    solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
end

"""
    demo_budget_controllers(; Q, T, delays)

Compare four budget spending controllers under information delay:
  - Naive (baseline)
  - Corrected denominator
  - Smith predictor
  - PID

One subplot per delay. Saves `budget_controllers.png`.
"""
function demo_budget_controllers(;
    Q = 100.0, T = 10.0,
    delays = [0.1, 0.3, 0.5] .* T,
    Kp = 1.0, Ki = 0.5, Kd = 0.1
)
    ts = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    plts = map(delays) do τ
        τ_pct = round(100.0 * τ / T; digits = 1)

        p = plot(ts, ideal_spent;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent",
            title = "τ = $(τ_pct)%·T",
            legend = :topleft, ylims = (0, :auto))
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)

        function add!(solver, label, color; idxs = 1, transform = s -> Q .- s)
            sol = solver(; Q, T, τ)
            ys = transform(sol.(ts; idxs = idxs))
            spent_end = round(ys[end]; digits = 1)
            plot!(p, ts, ys; label = "$label ($(spent_end))", linewidth = 2, color = color)
        end

        # scale gains by τ/T so PID stays stable across delay magnitudes
        τ_scaled_Kp = Kp / (1 + τ / T)
        τ_scaled_Ki = Ki / (1 + τ / T)^2
        τ_scaled_Kd = Kd / (1 + τ / T)
        pid_solver(; Q, T, τ) = solve_budget_pid(; Q, T, τ,
            Kp = τ_scaled_Kp, Ki = τ_scaled_Ki, Kd = τ_scaled_Kd)

        add!(solve_budget_delay,           "Naive",       :red)
        add!(solve_budget_corrected_denom, "Corr. denom", :blue)
        add!(solve_budget_smith,           "Smith",       :green;
             idxs = 2, transform = identity)
        add!(pid_solver,                   "PID",         :purple)
        p
    end

    fig = plot(plts...; layout = (length(delays), 1),
               size = (800, 360 * length(delays)),
               plot_title = "Budget spending controllers (Q=$Q, T=$T)")
    savefig(fig, "plots/budget_controllers.png")
    println("Plot saved to plots/budget_controllers.png")
    fig
end

"""
    demo_budget_controllers_noise(; Q, T, delays, noise_levels, n_samples, Kp, Ki, Kd)

Compare all four controllers (Naive, Corrected denom, Smith, PID) under
random τ noise.  Layout: one row per delay, one column per controller.
Each cell shows mean ± 1σ bands across `n_samples` trajectories with τ
drawn from Uniform(τ·(1-noise), τ·(1+noise)) for each noise level.

Saves `budget_controllers_noise.png`.
"""
function demo_budget_controllers_noise(;
    Q = 100.0, T = 10.0,
    delays      = [0.1, 0.3] .* T,
    noise_levels = [0.0, 0.10, 0.30],
    n_samples   = 40,
    Kp = 1.0, Ki = 0.5, Kd = 0.1
)
    ts = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    strategies = [
        ("Naive",           solve_budget_delay,           u -> Q .- u.(ts; idxs = 1)),
        ("Corr. denom",     solve_budget_corrected_denom, u -> Q .- u.(ts; idxs = 1)),
        ("Smith",           solve_budget_smith,           u -> u.(ts; idxs = 2)),
        ("PID",             (;Q,T,τ) -> solve_budget_pid(;Q,T,τ,
                                Kp=Kp/(1+τ/T), Ki=Ki/(1+τ/T)^2, Kd=Kd/(1+τ/T)),
                                                          u -> Q .- u.(ts; idxs = 1)),
    ]

    plts = []
    for τ in delays
        τ_pct = round(100.0 * τ / T; digits = 1)
        for (name, solver, extractor) in strategies
            p = plot(ts, ideal_spent;
                label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
                xlabel = "time", ylabel = "spent",
                title = "$name  (τ=$(τ_pct)%·T)",
                legend = :topleft, ylims = (0, :auto))
            hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)

            for noise in noise_levels
                cols = Vector{Vector{Float64}}()
                while length(cols) < max(1, n_samples)
                    τ_i = noise > 0 ? τ * (1 + noise * (2*rand() - 1)) : τ
                    sol = solver(; Q, T, τ = τ_i)
                    sol.retcode == ReturnCode.Success || continue
                    push!(cols, extractor(sol))
                end
                utils = hcat(cols...)
                μ = vec(mean(utils; dims = 2))
                σ = noise > 0 ? vec(std(utils; dims = 2)) : zeros(length(ts))
                spent_end = round(μ[end]; digits = 1)
                lbl = noise == 0.0 ? "no noise ($(spent_end))" :
                                     "±$(round(Int, noise*100))% ($(spent_end))"
                plot!(p, ts, μ; ribbon = σ, fillalpha = 0.15, linewidth = 2, label = lbl)
            end
            push!(plts, p)
        end
    end

    nrows = length(delays)
    fig = plot(plts...; layout = (nrows, length(strategies)),
               size = (380 * length(strategies), 400 * nrows),
               plot_title = "Controllers under random τ noise  (Q=$Q, T=$T)")
    savefig(fig, "plots/budget_controllers_noise.png")
    println("Plot saved to plots/budget_controllers_noise.png")
    fig
end

"""
    demo_budget_delay(; Q, T, delays)

Plot budget trajectories for several delay values alongside the ideal (τ=0)
pacing curve.  Saves `budget_delay.png`.

The ideal spend is B_ideal(t) = Q·(1 - t/T) (straight line to zero).
Each trajectory with τ > 0 shows how the controller over-estimates its
remaining budget and must accelerate spending as the deadline approaches.
"""
function demo_budget_delay(;
    Q = 100.0, T = 10.0,
    delays = [0.05, 0.1, 0.5] .* T,   # 5%, 10%, 50% of T
    tau_noise = 0.0,
    n_samples = 30
)
    ts = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    # one subplot per delay
    plts = map(delays) do τ
        τ_pct = round(100.0 * τ / T; digits = 1)

        p = plot(ts, ideal_spent;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent",
            title = "τ = $(τ_pct)% · T",
            legend = :topleft, ylims = (0, :auto))
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)

        function add_strategy!(p, solver, label, color)
            cols = Vector{Vector{Float64}}()
            while length(cols) < max(1, n_samples)
                τ_i = tau_noise > 0 ? τ * (1 + tau_noise * (2*rand() - 1)) : τ
                sol = solver(; Q, T, τ = τ_i)
                sol.retcode == ReturnCode.Success || continue
                push!(cols, Q .- sol.(ts; idxs = 1))
            end
            utils = hcat(cols...)
            μ = vec(mean(utils; dims = 2))
            σ = tau_noise > 0 ? vec(std(utils; dims = 2)) : zeros(length(ts))
            spent_end = round(μ[end]; digits = 1)
            plot!(p, ts, μ; ribbon = σ, fillalpha = 0.2, linewidth = 2,
                  color = color, label = "$label (spent: $spent_end)")
        end

        add_strategy!(p, solve_budget_delay,            "naive",     :red)
        add_strategy!(p, solve_budget_corrected_denom,  "corr. denom", :blue)
        add_strategy!(p, solve_budget_smith,            "Smith",     :green)
        p
    end

    title = "Budget spending strategies under information delay\n(Q=$Q, T=$T" *
            (tau_noise > 0 ? ", τ noise=$(round(Int, tau_noise*100))%" : "") * ")"
    fig = plot(plts...; layout = (length(delays), 1),
               size = (800, 350 * length(delays)),
               plot_title = title)
    savefig(fig, "plots/budget_delay.png")
    println("Plot saved to plots/budget_delay.png")
    fig
end

"""
    demo_budget_delay_with_noise(; Q, T, delays, noise_levels, n_samples)

Show how τ noise affects each of the three spending strategies.

Layout: one row per delay value, three columns (one per strategy).
Within each subplot several noise levels are overlaid as mean ± 1σ bands,
so the reader can see both the bias (mean deviation from ideal) and the
spread (sensitivity to τ uncertainty).

Saves `budget_delay_noise.png`.
"""
function demo_budget_delay_with_noise(;
    Q = 100.0, T = 10.0,
    delays = [0.1, 0.5] .* T,             # 10% and 50% of T
    noise_levels = [0.0, 0.05, 0.15, 0.30],
    n_samples = 40
)
    ts = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    strategies = [
        ("Naive",           solve_budget_delay,           u -> Q .- u.(ts; idxs = 1)),
        ("Corrected denom", solve_budget_corrected_denom, u -> Q .- u.(ts; idxs = 1)),
        ("Smith predictor", solve_budget_smith,           u -> u.(ts; idxs = 2)),
    ]

    plts = []
    for τ in delays
        τ_pct = round(100.0 * τ / T; digits = 1)
        for (col, (name, solver, extractor)) in enumerate(strategies)
            p = plot(ts, ideal_spent;
                label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
                xlabel = "time", ylabel = "spent",
                title = "$name  (τ=$(τ_pct)%·T)",
                legend = :topleft, ylims = (0, :auto))
            hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)

            for noise in noise_levels
                cols = Vector{Vector{Float64}}()
                while length(cols) < max(1, n_samples)
                    τ_i = noise > 0 ? τ * (1 + noise * (2*rand() - 1)) : τ
                    sol = solver(; Q, T, τ = τ_i)
                    sol.retcode == ReturnCode.Success || continue
                    push!(cols, extractor(sol))
                end
                utils = hcat(cols...)
                μ = vec(mean(utils; dims = 2))
                σ = noise > 0 ? vec(std(utils; dims = 2)) : zeros(length(ts))
                spent_end = round(μ[end]; digits = 1)
                lbl = noise == 0.0 ? "no noise ($(spent_end))" :
                                     "±$(round(Int, noise*100))% ($(spent_end))"
                plot!(p, ts, μ; ribbon = σ, fillalpha = 0.15, linewidth = 2, label = lbl)
            end
            push!(plts, p)
        end
    end

    nrows = length(delays)
    fig = plot(plts...; layout = (nrows, 3), size = (1100, 380 * nrows),
               plot_title = "Effect of τ noise on spending strategies  (Q=$Q, T=$T)")
    savefig(fig, "plots/budget_delay_noise.png")
    println("Plot saved to plots/budget_delay_noise.png")
    fig
end
