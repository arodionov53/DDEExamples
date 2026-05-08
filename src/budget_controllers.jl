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
    solve_budget_smith_adaptive(; Q, T, τ, α)

Adaptive Smith predictor: estimates the grant fulfilment rate α̂(t) from
delayed observations and uses it to scale the internal-model correction.

The standard Smith predictor computes:

    B̂(t) = B(t-τ) − (S(t) − S(t-τ))

and fails when the true fulfilment rate α < 1, because S(t) over-counts real
spend.  The adaptive version estimates α̂ from the ratio of observed drain to
predicted drain over the window [t-2τ, t-τ]:

    predicted drain = S(t-τ) − S(t-2τ)     (grants issued in that window)
    observed drain  = B(t-2τ) − B(t-τ)     (balance change seen via observation)
    α̂(t)           = observed / predicted  (clamped to [0.01, 2] for stability)

Then the corrected prediction is:

    B̂(t) = B(t-τ) − α̂(t) · (S(t) − S(t-τ))

When α̂ = 1 this reduces to the standard Smith predictor.  When α̂ < 1 the
controller issues grants at a higher rate to compensate.

State: u = [B(t), S(t)] — same as solve_budget_smith.
Requires three constant lags (τ, 2τ, 3τ) to estimate α̂.

`α` sets the true fulfilment rate of the simulated system (as in
solve_budget_smith_mismatch).  The controller does *not* know α — it
estimates α̂ from observations.
"""
function solve_budget_smith_adaptive(;
    Q = 100.0, T = 10.0, τ = 1.5, α = 1.0,
    tspan = (0.0, T - 1e-3)
)
    function smith_adaptive!(du, u, h, p, t)
        Q, T, τ, α = p
        dt = τ > 0 ? τ : 1e-6

        B_now   = h(p, t -   τ)[1];  S_now   = h(p, t -   τ)[2]
        B_prev  = h(p, t - 2τ)[1];  S_prev  = h(p, t - 2τ)[2]
        B_prev2 = h(p, t - 3τ)[1];  S_prev2 = h(p, t - 3τ)[2]

        # Estimate α̂ from the window [t-3τ, t-2τ]
        predicted_drain = S_prev - S_prev2          # grants issued
        observed_drain  = B_prev2 - B_prev          # actual balance drop
        α_hat = predicted_drain > 1e-6 ?
                    clamp(observed_drain / predicted_drain, 0.01, 2.0) : 1.0

        # Adaptive Smith prediction
        B_hat = B_now - α_hat * (u[2] - S_now)
        B_hat = max(B_hat, 0.0)

        rate  = B_hat / (T - t)
        du[1] = -α * rate       # actual spend
        du[2] =      rate       # controller's tally (grants issued)
    end

    h(p, t) = [Q, 0.0]
    p = (Q, T, τ, α)
    prob = DDEProblem(smith_adaptive!, [Q, 0.0], h, tspan, p;
                      constant_lags = [τ, 2τ, 3τ])
    solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
end

"""
    demo_smith_mismatch(; Q, T, τ, alphas)

Smith predictor under model mismatch: only a fraction α ∈ (0,1] of each
granted unit is actually consumed.  The controller's internal tally S(t)
counts every grant as fully consumed, but real spend is α × grant rate.

State: u = [B(t), S(t)] where S(t) is the controller's *believed* cumulative
spend (= α⁻¹ times actual spend).  The true balance evolves at rate α×rate,
while the controller predicts using S(t) at face value, causing it to believe
the budget is draining faster than it really is and therefore underspend.

α = 1.0  → perfect model, identical to solve_budget_smith (B(T) = 0).
α < 1.0  → model over-counts spend; controller underspends by ~(1-α).
"""
function solve_budget_smith_mismatch(;
    Q = 100.0, T = 10.0, τ = 1.5, α = 0.6,
    tspan = (0.0, T - 1e-3)
)
    function smith_mismatch!(du, u, h, p, t)
        Q, T, τ, α = p
        B_delayed = h(p, t - τ)[1]
        S_delayed = h(p, t - τ)[2]
        # Controller predicts balance using its (inflated) tally
        B_hat = B_delayed - (u[2] - S_delayed)
        rate  = max(B_hat, 0.0) / (T - t)   # grant rate computed by controller
        du[1] = -α * rate    # actual spend: only α fraction consumed
        du[2] =      rate    # controller's tally: counts full grant
    end

    h(p, t) = [Q, 0.0]
    p = (Q, T, τ, α)
    prob = DDEProblem(smith_mismatch!, [Q, 0.0], h, tspan, p; constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

"""
    demo_smith_mismatch(; Q, T, τ, alphas)

Show how the Smith predictor degrades when its internal model over-counts
actual spend (grant fulfilment rate α < 1).

One subplot per α value.  Each compares the Smith predictor (perfect model),
the mismatched Smith, and the corrected-denominator fallback.

Saves `plots/smith_mismatch.png`.
"""
function demo_smith_mismatch(;
    Q      = 100.0,
    T      = 10.0,
    τ      = 1.5,
    alphas = [0.9, 0.7, 0.5]
)
    ts          = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    plts = map(alphas) do α
        sol_perfect  = solve_budget_smith(; Q, T, τ)
        sol_mismatch = solve_budget_smith_mismatch(; Q, T, τ, α)
        sol_adaptive = solve_budget_smith_adaptive(; Q, T, τ, α)
        sol_cd       = solve_budget_corrected_denom(; Q, T, τ)

        spent_perfect  = sol_perfect.(ts;  idxs = 2)
        spent_mismatch = Q .- sol_mismatch.(ts; idxs = 1)
        spent_adaptive = Q .- sol_adaptive.(ts; idxs = 1)
        spent_cd       = Q .- sol_cd.(ts;  idxs = 1)

        p = plot(ts, ideal_spent;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent",
            title  = "α = $α  (τ = $(round(Int, 100τ/T))%·T)",
            legend = :topleft)
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)
        plot!(p, ts, spent_perfect;
            label = "Smith α=1 ($(round(spent_perfect[end]; digits=1)))",
            linewidth = 2, color = :green)
        plot!(p, ts, spent_mismatch;
            label = "Smith naive ($(round(spent_mismatch[end]; digits=1)))",
            linewidth = 2, color = :red)
        plot!(p, ts, spent_adaptive;
            label = "Smith adaptive ($(round(spent_adaptive[end]; digits=1)))",
            linewidth = 2, color = :purple)
        plot!(p, ts, spent_cd;
            label = "Corr. denom ($(round(spent_cd[end]; digits=1)))",
            linewidth = 2, color = :blue)
        p
    end

    fig = plot(plts...; layout = (length(alphas), 1),
               size = (800, 360 * length(alphas)),
               plot_title = "Smith predictor: naive vs. adaptive vs. corrected-denom  (Q=$Q, T=$T, τ=$τ)")
    savefig(fig, "plots/smith_mismatch.png")
    println("Plot saved to plots/smith_mismatch.png")
    fig
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
