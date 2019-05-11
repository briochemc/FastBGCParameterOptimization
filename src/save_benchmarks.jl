# This script it to run some benchmarks on 

# Load the packages needed
using BenchmarkTools, Random, JLD2

# List of functions to be benchmarked
list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_functions = [Symbol(string(m) * "_" * string(∇ᵏf̂)) for m in list_methods for ∇ᵏf̂ in list_∇ᵏf̂]

# Create a benchmark suite
suite = BenchmarkGroup()
myseed = rand(1:10000)
for fun in list_functions
    suite[string(fun)] = BenchmarkGroup() # Create a BenchmarkGroup for each function
    suite[string(fun)] = @benchmarkable $fun(λ₀ .+ randn(m)/10) seconds=18000 samples=1 evals=20 setup=(Random.seed!(myseed)) # Add benchmarkable objects in each group
    @eval $fun(λ₀ .* (1 .+ randn(m)/10)) # Run each function once for precompiling before benchmark
end

# Run the benchmarks
results = run(suite)

# Save the results in a Julia data file
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@save jld_file results list_functions

