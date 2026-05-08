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
