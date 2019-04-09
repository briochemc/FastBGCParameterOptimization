# This script is to run the optimization and print out the times as it goes in order to plot it after


# List of functions to be benchmarked
list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]

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

for (i, myrun) in enumerate(myruns), m in list_methods

    # objective, gradient, and Hessian symbols
    m_f̂t = Symbol(string(m) * "_f̂t")
    m_∇f̂t = Symbol(string(m) * "_∇f̂t")
    m_∇²f̂t = Symbol(string(m) * "_∇²f̂t")

    # Ensure the starting point is the same for everyone
    if m == :F1
        F1buf = F1.initialize_buffer(x₀, p₀)
    else
        mSol = Symbol(string(m) * "Sol")
        @eval $mSol = $m.Solution(copy(x₀))
    end

    println("\n┌────────────────────────")
    println("│ Optimizing using method " * string(m) * ", for " * myrun * "\n│")
    println("│    ", time())
    eval( :( results = optimize($m_f̂t, $m_∇f̂t, $m_∇²f̂t, $λ₀, NewtonTrustRegion(), $opt)) )
    println("└────────────────────────")
    # Save output
    jld_file = joinpath(path_to_package_root, "data", "Optim_callback_data" * string(m) * myruns2[i] * str_out * ".jld2")
    @save jld_file results list_methods
end



