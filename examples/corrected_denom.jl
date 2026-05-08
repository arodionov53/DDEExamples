# Corrected Denominator Budget Controller
#
# Problem: spend budget Q uniformly over horizon T, but the balance
# observation is delayed by τ.  The naive controller uses:
#
#   dB/dt = -B(t-τ) / (T - t)
#
# and overspends because B(t-τ) > B(t) always.
#
# Fix: replace (T-t) with (T-t+τ) in the denominator.
#
#   dB/dt = -B(t-τ) / (T - t + τ)
#
# Derivation: the observed balance B(t-τ) was correct τ time units ago.
# At that moment the remaining time was (T - (t-τ)) = (T - t + τ).
# Dividing by that larger number gives the rate that would have exhausted
# the budget *from the observation point*, which exactly compensates for
# the stale reading.
#
# Analytical check (τ=0): reduces to dB/dt = -B(t)/(T-t), exact solution
# B(t) = Q(1-t/T) — perfect linear drawdown.

using DelayDiffEq
using DifferentialEquations: Tsit5, ODEProblem
using Plots

# ── Parameters ────────────────────────────────────────────────────────────────

Q = 100.0   # total budget
T = 10.0    # horizon
delays = [0.0, 0.5, 1.5, 3.0]   # τ values to compare

ts = range(0.0, T - 1e-3; length = 500)
ideal = Q .* (1 .- ts ./ T)   # exact analytical solution at τ=0

# ── Solvers ───────────────────────────────────────────────────────────────────

function solve_naive(τ)
    if τ == 0.0
        prob = ODEProblem((u, p, t) -> [-u[1] / (T - t)], [Q], (0.0, T - 1e-3))
        return solve(prob, Tsit5())
    end
    function f!(du, u, h, p, t)
        du[1] = -h(p, t - τ)[1] / (T - t)
    end
    prob = DDEProblem(f!, [Q], (p, t) -> [Q], (0.0, T - 1e-3), ();
                      constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

function solve_corrected(τ)
    if τ == 0.0
        # τ=0: denominator becomes (T-t+0) = (T-t), same as naive → exact
        prob = ODEProblem((u, p, t) -> [-u[1] / (T - t)], [Q], (0.0, T - 1e-3))
        return solve(prob, Tsit5())
    end
    function f!(du, u, h, p, t)
        du[1] = -h(p, t - τ)[1] / (T - t + τ)   # ← corrected denominator
    end
    prob = DDEProblem(f!, [Q], (p, t) -> [Q], (0.0, T - 1e-3), ();
                      constant_lags = [τ])
    solve(prob, MethodOfSteps(Tsit5()))
end

# ── Plot: remaining balance B(t) ─────────────────────────────────────────────

colors_naive     = [:grey, :red,    :darkorange, :darkred]
colors_corrected = [:black, :blue,  :teal,       :purple]

p1 = plot(; xlabel = "time", ylabel = "remaining budget B(t)",
            title  = "Remaining balance",
            legend = :topright)
plot!(p1, ts, ideal; linewidth = 2, linestyle = :dash,
      color = :black, label = "ideal (τ=0)")

for (i, τ) in enumerate(delays)
    τ_pct = round(Int, 100 * τ / T)
    label = τ == 0.0 ? "" : "naive τ=$(τ_pct)%·T"
    sol = solve_naive(τ)
    τ == 0.0 && continue   # τ=0 naive = ideal, skip to avoid duplicate
    plot!(p1, ts, sol.(ts; idxs=1);
          color = colors_naive[i], linewidth = 1.5,
          linestyle = :dot, label = label)
end

for (i, τ) in enumerate(delays)
    τ_pct = round(Int, 100 * τ / T)
    label = τ == 0.0 ? "corr. τ=0 (ideal)" : "corr. τ=$(τ_pct)%·T"
    sol = solve_corrected(τ)
    plot!(p1, ts, sol.(ts; idxs=1);
          color = colors_corrected[i], linewidth = 2, label = label)
end

hline!(p1, [0.0]; linewidth = 1, linestyle = :dash, color = :grey, label = "")

# ── Plot: cumulative spent Q - B(t) ──────────────────────────────────────────

p2 = plot(; xlabel = "time", ylabel = "cumulative spent",
            title  = "Cumulative spend vs ideal",
            legend = :topleft)
plot!(p2, ts, Q .- ideal; linewidth = 2, linestyle = :dash,
      color = :black, label = "ideal")
hline!(p2, [Q]; linewidth = 1, linestyle = :dot, color = :grey, label = "cap Q")

for (i, τ) in enumerate(delays)
    τ == 0.0 && continue
    τ_pct = round(Int, 100 * τ / T)
    sol_n = solve_naive(τ)
    sol_c = solve_corrected(τ)
    spent_n = Q .- sol_n.(ts; idxs=1)
    spent_c = Q .- sol_c.(ts; idxs=1)
    final_n = round(spent_n[end]; digits=1)
    final_c = round(spent_c[end]; digits=1)
    plot!(p2, ts, spent_n; color = colors_naive[i],     linewidth = 1.5,
          linestyle = :dot, label = "naive τ=$(τ_pct)%·T ($final_n)")
    plot!(p2, ts, spent_c; color = colors_corrected[i], linewidth = 2,
          label = "corr. τ=$(τ_pct)%·T ($final_c)")
end

# ── Print summary table ───────────────────────────────────────────────────────

println("\nFinal spend as % of Q  (Q=$Q, T=$T)")
println("-"^48)
println("   τ       naive      corrected")
println("-"^48)
for τ in delays
    τ_pct = round(Int, 100 * τ / T)
    if τ == 0.0
        println("τ =  0%·T   100.0%      100.0%  (exact)")
    else
        sol_n = solve_naive(τ)
        sol_c = solve_corrected(τ)
        final_n = round(100 * (Q - sol_n(T - 1e-3; idxs=1)) / Q; digits=1)
        final_c = round(100 * (Q - sol_c(T - 1e-3; idxs=1)) / Q; digits=1)
        println("τ = $(lpad(τ_pct,2))%·T   $(lpad(final_n,6))%      $(lpad(final_c,6))%")
    end
end
println("-"^48)

# ── Save figure ───────────────────────────────────────────────────────────────

fig = plot(p1, p2; layout = (1, 2), size = (1000, 420),
           plot_title = "Corrected denominator controller  (Q=$Q, T=$T)")
savefig(fig, "corrected_denom.png")
println("\nPlot saved to corrected_denom.png")
