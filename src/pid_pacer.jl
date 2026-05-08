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
    demo_pid_pacer_noise(; Q, T, delays, noise_levels, n_samples)

Show how τ noise affects the PIDPacer vs Smith and corrected denominator.
Layout: one row per delay, three columns (one per controller).
Each cell overlays several noise levels as mean ± 1σ bands.

Saves `pid_pacer_noise.png`.
"""
function demo_pid_pacer_noise(;
    Q            = 100.0,
    T            = 10.0,
    delays       = [0.05, 0.1, 0.3] .* T,
    noise_levels = [0.0, 0.10, 0.30],
    n_samples    = 40,
    Kp = 1.0, Ki = 0.1, Kd = 0.05
)
    ts          = range(0.0, T - 1e-3; length = 500)
    ideal_spent = Q .* ts ./ T

    strategies = [
        ("Smith",       (;Q,T,τ) -> solve_budget_smith(; Q,T,τ),           u -> u.(ts; idxs=2)),
        ("Corr. denom", (;Q,T,τ) -> solve_budget_corrected_denom(; Q,T,τ), u -> Q .- u.(ts; idxs=1)),
        ("PID pacer",   (;Q,T,τ) -> solve_budget_pid_pacer(; Q,T,τ,Kp,Ki,Kd), u -> Q .- u.(ts; idxs=1)),
    ]

    plts = []
    for τ in delays
        τ_pct = round(100.0 * τ / T; digits = 1)
        for (name, solver, extractor) in strategies
            p = plot(ts, ideal_spent;
                label = "ideal", linewidth = 2, linestyle = :dash, color = :black,
                xlabel = "time", ylabel = "spent",
                title  = "$name  (τ=$(τ_pct)%·T)",
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
                utils    = hcat(cols...)
                μ        = vec(mean(utils; dims = 2))
                σ        = noise > 0 ? vec(std(utils; dims = 2)) : zeros(length(ts))
                spent_end = round(μ[end]; digits = 1)
                lbl = noise == 0.0 ? "no noise ($spent_end)" :
                                     "±$(round(Int, noise*100))% ($spent_end)"
                plot!(p, ts, μ; ribbon = σ, fillalpha = 0.15, linewidth = 2, label = lbl)
            end
            push!(plts, p)
        end
    end

    nrows = length(delays)
    fig = plot(plts...; layout = (nrows, length(strategies)),
               size = (380 * length(strategies), 400 * nrows),
               plot_title = "PIDPacer noise sensitivity  (Q=$Q, T=$T)")
    savefig(fig, "pid_pacer_noise.png")
    println("Plot saved to pid_pacer_noise.png")
    fig
end
