# This script it to run some benchmarks on 

# Load the packages needed
using BenchmarkTools, Random, JLD2

# Create a benchmark suite
suite = BenchmarkGroup()

# List of functions to be benchmarked
list_functions = [
    :q!
    :Dq!
    :D2q!
    :CSDDq!
    :FDDq!
    :ADDq!
]

# Create a BenchmarkGroup for each function
for fun in list_functions
    suite[string(fun)] = BenchmarkGroup()
end


myseed = rand(1:10000)

# Add some benchmarkable objects in each group
for fun in list_functions # for each function
    #suite[string(fun)] = @benchmarkable $fun(λᵣ) setup=(λᵣ = copy(λ₀ .* (1 .+ randn(npopt)/10)))seconds=18000
    #suite[string(fun)] = @benchmarkable $fun(λᵣ) seconds=18000 samples=2 setup=(λᵣ = copy(λ₀ .* (1 .+ randn(npopt)/10)))
    suite[string(fun)] = @benchmarkable $fun(λ₀ .+ randn(npopt)/10) seconds=18000 samples=1 evals=20 setup=(Random.seed!(myseed))
end

# Run each function once for precompiling before benchmark
for fun in list_functions # for each function
    @eval $fun(λ₀ .* (1 .+ randn(npopt)/10))
end

# Run the benchmarks
results = run(suite)

# Save the results in a Julia data file
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@save jld_file results list_functions

