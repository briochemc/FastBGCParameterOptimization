# This script is to run the optimization and print out the times as it goes in order to plot it after


# List of functions to be benchmarked
list_methods = [
    ("OF1"   , :f̂!,   :OF1_∇f̂!,   :OF1_∇²f̂!)
    ("F0"    , :f̂!,       :∇f̂!,    :F0_∇²f̂!)
    ("F1"    , :f̂!,       :∇f̂!,    :F1_∇²f̂!)
    ("DUAL"  , :f̂!,       :∇f̂!,  :DUAL_∇²f̂!)
    ("CSD"   , :f̂!,       :∇f̂!,   :CSD_∇²f̂!)
    ("FD1"   , :f̂!,       :∇f̂!,   :FD1_∇²f̂!)
    ("HYPER" , :f̂!, :HYPER_∇f̂!, :HYPER_∇²f̂!)
    ("FD2"   , :f̂!,   :FD2_∇f̂!,   :FD2_∇²f̂!)
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
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false, callback=print_full_state)

myruns = ["Compiling run", "Precompiled run"]
myruns2 = ["_compiling_run", ""]

using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)

for (i, myrun) in enumerate(myruns)
    for (method_name, f̂, ∇f̂, ∇²f̂) in list_methods
        println("\n┌────────────────────────")
        println("│ Optimizing using method " * method_name * ", for " * myrun * "\n│")
        init.x, init.p = 1x₀, 3p₀
        J.fac, J.p = factorize(∇ₓF(x₀, 3p₀)), 3p₀
        println("│    ", time())
        eval( :( results = optimize($f̂, $∇f̂, $∇²f̂, $λ₀, NewtonTrustRegion(), $opt)) )
        println("└────────────────────────")
        # Save output
        jld_file = joinpath(path_to_package_root, "data", "Optim_callback_data" * method_name * myruns2[i] * str_out * ".jld2")
        @save jld_file results list_methods
    end
end



