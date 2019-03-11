# This script is to run the optimization and print out the times as it goes in order to plot it after


# List of functions to be benchmarked
list_methods = [
    ("FLASH"     , :q!,  :Dq!,   :D2q!)
    ("FiniteDiff", :q!,  :Dq!,  :FDDq!)
    ("Dual"      , :q!,  :Dq!,  :ADDq!)
    ("Complex"   , :q!,  :Dq!, :CSDDq!)
    ("HyperDual" , :q!, :ADq!,  :AD2q!)
    ("sHyperDual", :q!, :ADq!, :AD2Sq!)
    ("bHyperDual", :dq,  :no!,    :no!)
]

function print_time_and_q(λ)
    println("   ", time(), q!(λ))
    return false
end

function print_time()
    println("   ", time())
    return false
end


function print_full_state(x)
    x.iteration == 0 ? println("│    time                   iteration      q(λ)           |Dq(λ)|") : nothing
    @printf "│    %15.20g" time()
    print(x)
    return false
end

# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false, callback=print_full_state)

using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)

for (method_name, q, Dq, D2q) in list_methods
    println("\n┌────────────────────────")
    println("│ Optimizing using method " * method_name * ", and printing time and q :\n│")
    init.x, init.p = 1x₀, 3p₀
    J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
    println("│    ", time())
    if method_name == "bHyperDual"
        eval( :( results = optimize($q, $λ₀, NewtonTrustRegion(), $opt)) )
    else
        eval( :( results = optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)) )
    end
    println("└────────────────────────")
    # Save output
    jld_file = joinpath(path_to_package_root, "data", "Optim_callback_data" * method_name * str_out * ".jld2")
    @save jld_file results list_methods
end



