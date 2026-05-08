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
