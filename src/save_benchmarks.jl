# This script it to run some benchmarks on 

# Load the packages needed
using BenchmarkTools, Random, JLD2

# List of functions to be benchmarked
list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_functions = [Symbol(string(m) * "_" * string(∇ᵏf̂)) for m in list_methods for ∇ᵏf̂ in list_∇ᵏf̂]

using Distributions, LinearAlgebra

μ = AIBECS.LNμ(mean_pobs, variance_pobs/100)
σ² = AIBECS.LNσ²(mean_pobs, variance_pobs/100)
σ = sqrt.(σ²)
λ_dist = MvNormal(μ, Diagonal(σ))

# Create a benchmark suite
suite = BenchmarkGroup()
myseed = 1
for fun in list_functions
    suite[string(fun)] = BenchmarkGroup() # Create a BenchmarkGroup for each function
    suite[string(fun)] = @benchmarkable $fun(rand(λ_dist)) seconds=18000 samples=1 evals=20 setup=(Random.seed!(myseed)) # Add benchmarkable objects in each group
    @eval $fun(rand(λ_dist)) # Run each function once for precompiling before benchmark
end

# Run the benchmarks
results = run(suite)

# Save the results in a Julia data file
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@save jld_file results list_functions

