module FD2

using DiffEqBase

mutable struct Solution
    s
end

function update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    if ~(sol.s isa SteadyStateSolution) || p ≠ sol.s.prob.p
        sol.s isa SteadyStateSolution ? x = sol.s.u : x = sol.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p) # define problem
        sol.s = solve(prob, alg; options...) # update s (inner solver)
    end
end

function f̂(f, F, ∇ₓF, sol, p, alg; options...) # objective
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s = sol.s.u
    return f(s,p)
end

function ∇f̂(f, F, ∇ₓF, sol, p, alg; options...)
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s, m, h = sol.s.u, length(p), 1e-4
    G = zeros(1,m)       # preallocate
    for j in 1:m
        hⱼ = h * p[j]
        p₊ = p + hⱼ * e(j,m)
        p₋ = p - hⱼ * e(j,m)
        s₊ = solve(SteadyStateProblem(F, ∇ₓF, s, p₊), alg; options...).u
        s₋ = solve(SteadyStateProblem(F, ∇ₓF, s, p₋), alg; options...).u
        G[j] = (f(s₊,p₊) - f(s₋,p₋)) / 2hⱼ
    end
    return G
end

function ∇²f̂(f, F, ∇ₓF, sol, p, alg; options...) # Hessian
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s, m, h = sol.s.u, length(p), 1e-2
    H = zeros(m,m)       # preallocate
    for j in 1:m
        hⱼ = h * p[j]
        p₊ = p + hⱼ * e(j,m)
        p₋ = p - hⱼ * e(j,m)
        s₊ = solve(SteadyStateProblem(F, ∇ₓF, s, p₊), alg; options...).u
        s₋ = solve(SteadyStateProblem(F, ∇ₓF, s, p₋), alg; options...).u
        H[j,j] = (f(s₋,p₋) - 2f(s,p) + f(s₊,p₊)) / hⱼ^2
        for k = j+1:m # Off-diagonal terms
            hₖ = h * p[k]
            p₊₊ = p₊ + hₖ * e(k,m)
            p₊₋ = p₊ - hₖ * e(k,m)
            p₋₊ = p₋ + hₖ * e(k,m)
            p₋₋ = p₋ - hₖ * e(k,m)
            s₊₊ = solve(SteadyStateProblem(F, ∇ₓF, s₊, p₊₊), alg; options...).u
            s₊₋ = solve(SteadyStateProblem(F, ∇ₓF, s₊, p₊₋), alg; options...).u
            s₋₊ = solve(SteadyStateProblem(F, ∇ₓF, s₋, p₋₊), alg; options...).u
            s₋₋ = solve(SteadyStateProblem(F, ∇ₓF, s₋, p₋₋), alg; options...).u
            H[j,k] = (f(s₊₊,p₊₊) - f(s₋₊,p₋₊) - f(s₊₋,p₊₋) + f(s₋₋,p₋₋)) / (4 * hⱼ * hₖ)
            j ≠ k ? H[k,j] = H[j,k] : nothing
        end
    end
    return H
end

e(j, m) = [i == j for i in 1:m]      # j-th basis vector

end # module
