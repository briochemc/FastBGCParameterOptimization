# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false)

# Initialize TimerOutput
using TimerOutputs
const to = TimerOutput()

# Time all the functions
list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
for m in list_methods, ∇ᵏf̂ in list_∇ᵏf̂
    label = string(m) * "_" * string(∇ᵏf̂)
    m_∇ᵏf̂ = Symbol(label)
    m_∇ᵏf̂t = Symbol(label * "t")
    if ∇ᵏf̂ == :f̂
        @eval $m_∇ᵏf̂t(λ) = @timeit to $label $m_∇ᵏf̂(λ)
    else
        @eval $m_∇ᵏf̂t(s, λ) = @timeit to $label $m_∇ᵏf̂(s, λ)
    end
end

# path to root of directory for loading and saving
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)

# Two runs, 1st run to precompile
myruns = ["Compiling run", "Precompiled run"]
myruns2 = ["_compiling_run", ""]

# timers as a Dict to hold all timer data
timers = Dict()
function print_mytimer!(d::Dict, f::String, t::TimerOutput)
    d2 = Dict() # sub dictionary to fill with the timer data
    t2 = t.inner_timers["Trust Region"]
    push!(d2, "total_time" => copy(t2.accumulated_data.time))
    push!(d2, "total_mem" => copy(t2.accumulated_data.allocs))
    for f2 in keys(t2.inner_timers)
        t3 = t2.inner_timers[f2]
        d3 = Dict()
        push!(d3, "time" => copy(t3.accumulated_data.time))
        push!(d3, "ncalls" => copy(t3.accumulated_data.ncalls))
        push!(d3, "allocs" => copy(t3.accumulated_data.allocs))
        push!(d2, string(f2) => d3)
    end
    push!(d, f => d2)
end

# Run the optimizations
for (i, myrun) in enumerate(myruns), m in list_methods
    println("\n\n------------------------\n") # print what you are doing
    println(myrun * ": " * string(m) * " (with TimerOutputs)")
    println("\n------------------------\n\n")

    # Ensure the starting point is the same for everyone
    if m == :F1
        F1buf = F1.initialize_buffer(x₀, p₀)
    else
        mSol = Symbol(string(m) * "Sol")
        @eval $mSol = $m.Solution(copy(x₀))
    end

    # File path to print TimerOutputs to specific method and run
    file_name = joinpath(path_to_package_root, "data", "TimerOutputs_" * string(m) * myruns2[i] * str_out * ".txt")

    # objective, gradient, and Hessian symbols
    m_f̂t = Symbol(string(m) * "_f̂t")
    m_∇f̂t = Symbol(string(m) * "_∇f̂t")
    m_∇²f̂t = Symbol(string(m) * "_∇²f̂t")

    # Open file, run the optimization, and write to file
    io = open(file_name, "w")
        # Reset timer
        reset_timer!(to) 
        # Run the timed optimization
        @timeit to "Trust Region" eval(:(res = optimize($m_f̂t, $m_∇f̂t, $m_∇²f̂t, $λ₀, NewtonTrustRegion(), $opt)))
        # Print the TimerOutput
        print_timer(io, to)
        print(io, string(m) * "\n" * read(io, String))
    close(io)

    # Save into Dict
    print_mytimer!(timers, string(m) * myruns2[i], to)
end

# Save the results in a Julia data file
using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data", "TimerOutputs_data" * str_out * ".jld2")
@save jld_file timers list_methods



