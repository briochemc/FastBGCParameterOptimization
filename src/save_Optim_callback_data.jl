# This script is to run the optimization and print out the times as it goes in order to plot it after


# List of functions to be benchmarked
list_methods = [
    ("FD1"        , :f̂!,   :∇f̂!,  :FD∇f̂!)
    ("CSD"        , :f̂!,   :∇f̂!, :CSD∇f̂!)
    ("DUAL"       , :f̂!,   :∇f̂!,  :AD∇f̂!)
    ("FLASH"      , :f̂!,   :∇f̂!,   :∇²f̂!)
    ("HYPER"      , :f̂!,  :ADf̂!,  :AD²f̂!)
    ("HYPERSMART" , :f̂!, :HS∇f̂!, :HS∇²f̂!)
    ("FD2"        , :f̂!,  :FDf̂!,  :FD²f̂!)
]

function print_time_and_q(λ)
    println("   ", time(), f̂!(λ))
    return false
end

function print_time()
    println("   ", time())
    return false
end


function print_full_state(x)
    x.iteration == 0 ? println("│    time                   iteration      f̂(λ)           |∇f̂(λ)|") : nothing
    @printf "│    %20.16e" time()
    print(x)
    return false
end

# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false, callback=print_full_state, g_calls_limit=15)

using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)

for (method_name, f̂, ∇f̂, ∇²f̂) in list_methods
    println("\n┌────────────────────────")
    println("│ Optimizing using method " * method_name * ", and printing time:\n│")
    init.x, init.p = 1x₀, 3p₀
    J.fac, J.p = factorize(∇ₓF(x₀, 3p₀)), 3p₀
    println("│    ", time())
    eval( :( results = optimize($f̂, $∇f̂, $∇²f̂, $λ₀, NewtonTrustRegion(), $opt)) )
    println("└────────────────────────")
    # Save output
    jld_file = joinpath(path_to_package_root, "data", "Optim_callback_data" * method_name * str_out * ".jld2")
    @save jld_file results list_methods
end



