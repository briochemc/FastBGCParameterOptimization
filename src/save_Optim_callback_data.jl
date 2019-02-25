# This script is to run the optimization and print out the times as it goes in order to plot it after 


# List of functions to be benchmarked
list_methods = [
    (  "D2q", :q!,  :Dq!,   :D2q!)
    ( "FDDq", :q!,  :Dq!,  :FDDq!)
    ( "ADDq", :q!,  :Dq!,  :ADDq!)
    ("CSDDq", :q!,  :Dq!, :CSDDq!)
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
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false, x_tol = 1e-3, callback=print_full_state)

for (method_name, q, Dq, D2q) in list_methods
    println("\n┌────────────────────────")
    println("│ Optimizing using method " * method_name * ", and printing time and q :\n│")
    init.x, init.p = 1x₀, 3p₀
    J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
    eval( :( optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)) ) 
    println("└────────────────────────")
end


