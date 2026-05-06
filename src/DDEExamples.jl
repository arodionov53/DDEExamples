module DDEExamples

using DelayDiffEq
using DifferentialEquations: Tsit5, ODEProblem
using Plots
using Statistics: mean, std

export solve_mackey_glass, solve_logistic_dde, solve_two_delay, solve_random_delay,
       solve_mackey_glass_nodelay, solve_logistic_nodelay, solve_two_delay_nodelay,
       solve_budget_delay, solve_budget_corrected_denom, solve_budget_smith,
       solve_budget_pid, demo_budget_controllers, demo_budget_controllers_noise,
       solve_budget_pid_pacer, demo_pid_pacer,
       demo_budget_delay, demo_budget_delay_with_noise,
       demo, demo_zero_delay

"""
Mackey-Glass equation: du/dt = β * u(t-τ) / (1 + u(t-τ)^n) - γ * u(t)

A classic DDE from physiology (blood cell regulation) that exhibits chaotic
behavior for certain parameter values.
"""
function solve_mackey_glass(;
    β = 2.0, γ = 1.0, n = 9.65, τ = 2.0,
    tspan = (0.0, 100.0), u0 = 0.5
)
    function mackey_glass!(du, u, h, p, t)
        β, γ, n, τ = p
        hist = h(p, t - τ)
        du[1] = β * hist[1] / (1 + hist[1]^n) - γ * u[1]
    end

    h(p, t) = [0.5]
    p = (β, γ, n, τ)
    prob = DDEProblem(mackey_glass!, [u0], h, tspan, p; constant_lags = [τ])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    sol
end

"""
Logistic DDE: du/dt = r * u(t) * (1 - u(t-τ) / K)

Delayed logistic growth — the population's growth rate depends on the
population size τ time units in the past.
"""
function solve_logistic_dde(;
    r = 0.5, K = 10.0, τ = 5.0,
    tspan = (0.0, 80.0), u0 = 0.1
)
    function logistic!(du, u, h, p, t)
        r, K, τ = p
        hist = h(p, t - τ)
        du[1] = r * u[1] * (1 - hist[1] / K)
    end

    h(p, t) = [0.1]
    p = (r, K, τ)
    prob = DDEProblem(logistic!, [u0], h, tspan, p; constant_lags = [τ])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    sol
end

"""
Two-delay system: du/dt = -a * u(t-τ₁) - b * u(t-τ₂)

A system with two distinct delays, showing how multiple lags interact.
"""
function solve_two_delay(;
    a = 1.0, b = 0.5, τ₁ = 1.0, τ₂ = 3.0,
    tspan = (0.0, 40.0), u0 = 1.0
)
    function two_delay!(du, u, h, p, t)
        a, b, τ₁, τ₂ = p
        h1 = h(p, t - τ₁)
        h2 = h(p, t - τ₂)
        du[1] = -a * h1[1] - b * h2[1]
    end

    h(p, t) = [1.0]
    p = (a, b, τ₁, τ₂)
    prob = DDEProblem(two_delay!, [u0], h, tspan, p; constant_lags = [τ₁, τ₂])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    sol
end

"""
Logistic DDE with random delay: du/dt = r * u(t) * (1 - u(t-τ) / K), τ ~ Uniform[τ_min, τ_max]

Uses an EnsembleProblem to run many trajectories, each with a different delay
sampled uniformly at random.  This models uncertainty in the feedback lag —
for example, a biological maturation time that varies between individuals.
"""
function solve_random_delay(;
    r = 0.5, K = 10.0, τ_min = 2.0, τ_max = 8.0,
    tspan = (0.0, 80.0), u0 = 0.1, trajectories = 20
)
    function logistic!(du, u, h, p, t)
        r, K, τ = p
        hist = h(p, t - τ)
        du[1] = r * u[1] * (1 - hist[1] / K)
    end

    h(p, t) = [u0]
    τ_mid = (τ_min + τ_max) / 2
    p = (r, K, τ_mid)
    prob = DDEProblem(logistic!, [u0], h, tspan, p; constant_lags = [τ_mid])

    function prob_func(prob, ctx)
        τ_rand = τ_min + (τ_max - τ_min) * rand(ctx.rng)
        p_new = (prob.p[1], prob.p[2], τ_rand)
        remake(prob, p = p_new, constant_lags = [τ_rand])
    end

    ensemble = EnsembleProblem(prob; prob_func = prob_func)
    sim = solve(ensemble, MethodOfSteps(Tsit5()), EnsembleSerial(); trajectories = trajectories)
    sim
end

"""
Mackey-Glass with zero delay (ODE): du/dt = β * u / (1 + u^n) - γ * u

When τ = 0 the delayed term u(t-τ) becomes u(t) and the equation reduces to
an autonomous ODE that converges to a stable equilibrium.
"""
function solve_mackey_glass_nodelay(;
    β = 2.0, γ = 1.0, n = 9.65,
    tspan = (0.0, 100.0), u0 = 0.5
)
    function mg_ode!(du, u, p, t)
        β, γ, n = p
        du[1] = β * u[1] / (1 + u[1]^n) - γ * u[1]
    end

    p = (β, γ, n)
    prob = ODEProblem(mg_ode!, [u0], tspan, p)
    solve(prob, Tsit5())
end

"""
Logistic growth with zero delay (ODE): du/dt = r * u * (1 - u / K)

When τ = 0 the equation is the standard logistic ODE — the solution is a
smooth sigmoid that monotonically approaches the carrying capacity K.
"""
function solve_logistic_nodelay(;
    r = 0.5, K = 10.0,
    tspan = (0.0, 80.0), u0 = 0.1
)
    function logistic_ode!(du, u, p, t)
        r, K = p
        du[1] = r * u[1] * (1 - u[1] / K)
    end

    p = (r, K)
    prob = ODEProblem(logistic_ode!, [u0], tspan, p)
    solve(prob, Tsit5())
end

"""
Two-delay system with zero delay (ODE): du/dt = -(a + b) * u

When both delays are zero the equation reduces to exponential decay
u(t) = u₀ exp(-(a+b)t).
"""
function solve_two_delay_nodelay(;
    a = 1.0, b = 0.5,
    tspan = (0.0, 40.0), u0 = 1.0
)
    function td_ode!(du, u, p, t)
        a, b = p
        du[1] = -(a + b) * u[1]
    end

    p = (a, b)
    prob = ODEProblem(td_ode!, [u0], tspan, p)
    solve(prob, Tsit5())
end

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
    savefig(fig, "budget_controllers.png")
    println("Plot saved to budget_controllers.png")
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
    savefig(fig, "budget_controllers_noise.png")
    println("Plot saved to budget_controllers_noise.png")
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
    savefig(fig, "budget_delay.png")
    println("Plot saved to budget_delay.png")
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
    savefig(fig, "budget_delay_noise.png")
    println("Plot saved to budget_delay_noise.png")
    fig
end

"""
    solve_budget_pid_pacer(; Q, T, τ, Kp, Ki, Kd, max_integral, request_rate)

DDE model of the production PIDPacer from pid_pacer.go.

Faithfully reproduces the Go implementation:
  - Error signal is rate-based: e(t) = targetRate - observedRate(t-τ)
  - observedRate(t) = -dB/dt = spend rate inferred from delayed balance derivative
  - PID output mapped through sigmoid: prob(t) = 1 / (1 + exp(-PID))
  - Integral windup protection: integral clamped to ±max_integral/Ki
  - spend(t) = prob(t) * request_rate   (probabilistic grants)

State: u = [B(t), I(t)] where I is the clamped integral of the rate error.

Parameters match PIDConfig defaults: Kp=1.0, Ki=0.1, Kd=0.05.
`request_rate` is the average spend per second if all grants are accepted.
"""
function solve_budget_pid_pacer(;
    Q            = 100.0,
    T            = 10.0,
    τ            = 1.0,
    Kp           = 1.0,
    Ki           = 0.1,
    Kd           = 0.05,
    max_integral = 10.0 / 0.1,   # mirrors Go: 10.0 / Ki
    request_rate = Q / T,         # average grant size if prob=1
    tspan        = (0.0, T - 1e-3)
)
    target_rate = Q / T

    function pid_pacer!(du, u, h, p, t)
        Q, T, τ, Kp, Ki, Kd, max_integral, target_rate, request_rate = p

        # Observed spend rate from delayed balance: finite diff over τ
        B_now  = h(p, t - τ)[1]
        B_prev = h(p, t - 2τ)[1]
        dt_obs = τ > 0 ? τ : 1e-6
        observed_rate = (B_prev - B_now) / dt_obs   # positive = spending

        # Rate error (mirrors Go: targetRate - currentRate)
        e = target_rate - observed_rate

        # Integral with windup protection (mirrors Go clamping)
        I_new = clamp(u[2] + e, -max_integral, max_integral)

        # Derivative of error
        B_now2 = h(p, t - 3τ)[1]
        observed_rate_prev = (B_now - B_now2) / dt_obs
        e_prev = target_rate - observed_rate_prev
        de = (e - e_prev) / dt_obs

        # PID output → sigmoid probability (mirrors Go)
        pid_out = Kp * e + Ki * I_new + Kd * de
        prob    = 1.0 / (1.0 + exp(-pid_out))
        prob    = clamp(prob, 0.0, 1.0)

        # Spend rate = probability × request_rate
        du[1] = -prob * request_rate
        du[2] = e   # integral state (clamping applied above via I_new)
        # Note: we track unclamped I in u[2] for solver continuity,
        # but use clamped I_new in the PID calculation
    end

    h(p, t) = [Q, 0.0]
    p = (Q, T, τ, Kp, Ki, Kd, max_integral, target_rate, request_rate)
    prob = DDEProblem(pid_pacer!, [Q, 0.0], h, tspan, p;
                      constant_lags = [τ, 2τ, 3τ])
    solve(prob, MethodOfSteps(Tsit5()); dtmax = τ / 10)
end

"""
    demo_pid_pacer(; Q, T, delays, noise_levels, n_samples)

Compare the production PIDPacer DDE model against Smith predictor and
corrected denominator across several delay values and (optionally) noise.

One subplot per delay. Saves `pid_pacer_comparison.png`.
"""
function demo_pid_pacer(;
    Q            = 100.0,
    T            = 10.0,
    delays       = [0.05, 0.1, 0.3] .* T,
    tau_noise    = 0.0,
    n_samples    = 30,
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts           = range(0.0, T - 1e-3; length = 500)
    ideal_spent  = Q .* ts ./ T

    plts = map(delays) do τ
        τ_pct = round(100.0 * τ / T; digits = 1)

        p = plot(ts, ideal_spent;
            label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
            xlabel = "time", ylabel = "spent",
            title  = "τ = $(τ_pct)%·T",
            legend = :topleft, ylims = (0, :auto))
        hline!(p, [Q]; label = "cap", linewidth = 1, linestyle = :dot, color = :grey)

        function add!(label, color, solver; idxs=1, xform = s -> Q .- s)
            cols = Vector{Vector{Float64}}()
            while length(cols) < max(1, tau_noise > 0 ? n_samples : 1)
                τ_i = tau_noise > 0 ? τ * (1 + tau_noise*(2*rand()-1)) : τ
                sol = solver(τ_i)
                sol.retcode == ReturnCode.Success || continue
                push!(cols, xform(sol.(ts; idxs=idxs)))
            end
            utils = hcat(cols...)
            μ = vec(mean(utils; dims=2))
            σ = size(utils,2) > 1 ? vec(std(utils; dims=2)) : zeros(length(ts))
            spent_end = round(μ[end]; digits=1)
            plot!(p, ts, μ; ribbon=σ, fillalpha=0.15, linewidth=2,
                  color=color, label="$label ($spent_end)")
        end

        add!("Smith",       :green,  τ_i -> solve_budget_smith(; Q, T, τ=τ_i);
             idxs=2, xform=identity)
        add!("Corr. denom", :blue,   τ_i -> solve_budget_corrected_denom(; Q, T, τ=τ_i))
        add!("PID pacer",   :purple, τ_i -> solve_budget_pid_pacer(; Q, T, τ=τ_i, Kp, Ki, Kd))
        p
    end

    noise_str = tau_noise > 0 ? ", τ noise=±$(round(Int,tau_noise*100))%" : ""
    fig = plot(plts...; layout=(length(delays), 1),
               size=(800, 380*length(delays)),
               plot_title="PIDPacer vs Smith vs Corr.denom  (Q=$Q, T=$T$noise_str)")
    savefig(fig, "pid_pacer_comparison.png")
    println("Plot saved to pid_pacer_comparison.png")
    fig
end

"""
    demo_zero_delay()

Compare each DDE with its zero-delay (ODE) counterpart. Produces a plot saved
to `dde_vs_ode.png` showing how the delay fundamentally changes system behavior.
"""
function demo_zero_delay()
    sol_mg = solve_mackey_glass()
    sol_mg0 = solve_mackey_glass_nodelay()
    sol_lg = solve_logistic_dde()
    sol_lg0 = solve_logistic_nodelay()
    sol_td = solve_two_delay()
    sol_td0 = solve_two_delay_nodelay()

    p1 = plot(sol_mg; idxs = 1, label = "τ = 2 (DDE)", xlabel = "t", ylabel = "u(t)",
              title = "Mackey-Glass: delay vs no delay")
    plot!(p1, sol_mg0; idxs = 1, label = "τ = 0 (ODE)", linestyle = :dash, linewidth = 2)

    p2 = plot(sol_lg; idxs = 1, label = "τ = 5 (DDE)", xlabel = "t", ylabel = "u(t)",
              title = "Logistic Growth: delay vs no delay")
    plot!(p2, sol_lg0; idxs = 1, label = "τ = 0 (ODE)", linestyle = :dash, linewidth = 2)

    p3 = plot(sol_td; idxs = 1, label = "τ₁=1, τ₂=3 (DDE)", xlabel = "t", ylabel = "u(t)",
              title = "Two-Delay System: delay vs no delay")
    plot!(p3, sol_td0; idxs = 1, label = "τ = 0 (ODE)", linestyle = :dash, linewidth = 2)

    fig = plot(p1, p2, p3; layout = (3, 1), size = (800, 900))
    savefig(fig, "dde_vs_ode.png")
    println("Plot saved to dde_vs_ode.png")
    fig
end

"""
    demo()

Solve all four example DDEs and produce a combined plot saved to `dde_examples.png`.
"""
function demo()
    sol_mg = solve_mackey_glass()
    sol_lg = solve_logistic_dde()
    sol_td = solve_two_delay()
    sim_rd = solve_random_delay(; trajectories = 20)

    p1 = plot(sol_mg; idxs = 1, title = "Mackey-Glass (chaotic)", xlabel = "t", ylabel = "u(t)", legend = false)
    p2 = plot(sol_lg; idxs = 1, title = "Delayed Logistic Growth", xlabel = "t", ylabel = "u(t)", legend = false)
    p3 = plot(sol_td; idxs = 1, title = "Two-Delay System", xlabel = "t", ylabel = "u(t)", legend = false)

    p4 = plot(; title = "Random Delay Logistic (τ ∈ [2, 8])", xlabel = "t", ylabel = "u(t)", legend = false)
    for sol in sim_rd.u
        plot!(p4, sol.t, [u[1] for u in sol.u]; alpha = 0.4)
    end

    fig = plot(p1, p2, p3, p4; layout = (4, 1), size = (800, 1200))
    savefig(fig, "dde_examples.png")
    println("Plot saved to dde_examples.png")
    fig
end

end # module DDEExamples
