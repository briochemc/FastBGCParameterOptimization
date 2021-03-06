# This script is to run the optimization and save the states along the optimization

# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = true)

using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)

AF1_mem = F1.initialize_mem(x₀, p₀)

println("\n┌────────────────────────")
println("│ Optimizing using method AF1\n│")
println("│    ", time())
results = optimize(AF1_f̂, AF1_∇f̂, AF1_∇²f̂, λ₀, NewtonTrustRegion(), opt)
println("└────────────────────────")

# Save output
jld_file = joinpath(path_to_package_root, "data", "Optim_trace.jld2")

λopt = results.minimizer
popt = λ2p(λopt)
probopt = SteadyStateProblem(A_F, A_∇ₓF, x₀, popt)
xopt = solve(probopt, CTKAlg(), nrm=nrm).u

xs = zeros(length(x₀), length(results.trace))

for (i, trace_i) in enumerate(results.trace)
    λ = trace_i.metadata["x"]
    p = λ2p(λ)
    prob = SteadyStateProblem(A_F, A_∇ₓF, x₀, p)
    xs[:,i] .= solve(prob, CTKAlg(), nrm=nrm).u
end

@save jld_file results x₀ p₀ λ₀ xopt λopt popt μDIPobs xs



