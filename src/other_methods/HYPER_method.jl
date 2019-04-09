module HYPER

using LinearAlgebra, DiffEqBase
using HyperDualNumbers, HyperDualMatrixTools
using DualNumbers, DualMatrixTools

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
    s, m = sol.s.u, length(p)
    out = zeros(1,m)       # preallocate
    for j in 1:m
        pⱼ = p + ε * e(j,m)           # Dual p
        prob = SteadyStateProblem(F, ∇ₓF, s, pⱼ) # define problem
        sⱼ = solve(prob, alg; options...).u # update s (inner solver)
        out[j] = 𝔇(f(sⱼ, pⱼ)) # Dual pass
    end
    return out
end

function ∇²f̂(f, F, ∇ₓF, sol, p, alg; options...) # Hessian
    update_Solution!(F, ∇ₓF, sol, p, alg; options...)
    s, m = sol.s.u, length(p)
    out = zeros(m,m)       # preallocate
    for j in 1:m, k in j:m
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)     # Hyperdual p
        prob = SteadyStateProblem(F, ∇ₓF, s, pⱼₖ) # define problem
        sⱼₖ = solve(prob, alg; options...).u # update s (inner solver)
        out[j,k] = ℌ(f(sⱼₖ, pⱼₖ)) # HyperDual pass
        j ≠ k ? out[k,j] = out[j,k] : nothing
    end
    return out
end

# Helper functions
e(j, m) = [i == j for i in 1:m]      # j-th basis vector
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x)      # dual part

end # module
