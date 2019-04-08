module DUAL

using LinearAlgebra
using DualNumbers, DualMatrixTools, DiffEqBase

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

# Analytical Jacobian formula (transpose of Eq.(?))
∇f̂(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s, p) = ∇ₚf(s,p) - (∇ₚF(s,p)' * (∇ₓF(s,p)' \ ∇ₓf(s,p)'))'

function ∇f̂(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...)
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s = sol.s.u
    return ∇f̂(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, s, p)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...) # Hessian
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s, m = sol.s.u, length(p)
    out = zeros(m,m)       # preallocate
    for j in 1:m
        pⱼ = p + ε * e(j,m)           # Dual p
        prob = SteadyStateProblem(F, ∇ₓF, s, pⱼ) # define problem
        sⱼ = solve(prob, alg; options...).u # update s (inner solver)
        out[j,:] .= vec(𝔇(∇f̂(∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sⱼ, pⱼ))) # Dual of Eq.(?)
    end
    return out
end

# Helper functions
e(j, m) = [i == j for i in 1:m]      # j-th basis vector
𝔇(x) = DualNumbers.dualpart.(x)      # dual part

end # module
