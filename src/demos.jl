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
    savefig(fig, "plots/dde_vs_ode.png")
    println("Plot saved to plots/dde_vs_ode.png")
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
    savefig(fig, "plots/dde_examples.png")
    println("Plot saved to plots/dde_examples.png")
    fig
end
