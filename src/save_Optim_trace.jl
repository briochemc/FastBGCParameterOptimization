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
@save jld_file results x₀ p₀ λ₀



