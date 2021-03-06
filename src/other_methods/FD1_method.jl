module FD1

using LinearAlgebra, DiffEqBase

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

function objective(f, F, ∇ₓF, sol, p, alg; options...)
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s = sol.s.u
    return f(s,p)
end

# Analytical Jacobian formula (transpose of Eq.(14))
gradient(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s, p) = ∇ₚf(s,p) - (∇ₚF(s,p)' * (∇ₓF(s,p)' \ ∇ₓf(s,p)'))'

function gradient(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...)
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s = sol.s.u
    return gradient(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s, p)
end

function hessian(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...) 
    s, m, h = sol.s.u, length(p), 1e-4
    H = zeros(m,m)       # preallocate
    for j in 1:m
        hⱼ = h * p[j]
        p₊ = p + hⱼ * e(j,m)
        p₋ = p - hⱼ * e(j,m)
        s₊ = solve(SteadyStateProblem(F, ∇ₓF, s, p₊), alg; options...).u
        s₋ = solve(SteadyStateProblem(F, ∇ₓF, s, p₋), alg; options...).u
        H[j,:] .= vec((gradient(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s₊,p₊) - gradient(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s₋,p₋)) / 2hⱼ)
    end
    return H
end

e(j, m) = [i == j for i in 1:m]      # j-th basis vector

end # module
