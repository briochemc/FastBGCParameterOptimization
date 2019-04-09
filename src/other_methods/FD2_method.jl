module FD2

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

function ∇f̂(f, F, ∇ₓF, sol, p, alg; options...)
    fun = p -> [f̂(f, F, ∇ₓF, sol, p, alg; options...)]
    return Calculus.jacobian(fun, p, :central)
end

function ∇²f̂(f, F, ∇ₓF, sol, p, alg; options...) # Hessian
    fun = p -> f̂(f, F, ∇ₓF, sol, p, alg; options...)
    return Calculus.hessian(fun, p)
end

end # module
