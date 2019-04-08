module FD1

using LinearAlgebra
using Calculus, DiffEqBase

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

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...) 
    fun = p -> vec(∇f̂(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, sol, p, alg; options...))
    return Calculus.jacobian(fun, p, :central)
end

end # module
