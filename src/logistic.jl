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
