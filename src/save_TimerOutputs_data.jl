# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, extended_trace = true)
using TimerOutputs
const to = TimerOutput()

# timed version for all functions called in Optim
qt!(λ)        = @timeit to "f"   q!(λ)
Dqt!(s, λ)    = @timeit to "∇f"  Dq!(s, λ)
ADqt!(s, λ)   = @timeit to "∇f"  ADq!(s, λ)
FDqt!(s, λ)   = @timeit to "∇f"  FDq!(s, λ)
HSDqt!(s, λ)  = @timeit to "∇f"  HSDq!(s, λ)
D2qt!(s, λ)   = @timeit to "∇²f" D2q!(s, λ)
ADDqt!(s, λ)  = @timeit to "∇²f" ADDq!(s, λ)
AD2qt!(s, λ)  = @timeit to "∇²f" AD2q!(s, λ)
HSD2qt!(s, λ) = @timeit to "∇²f" HSD2q!(s, λ)
CSDDqt!(s, λ) = @timeit to "∇²f" CSDDq!(s, λ)
FDDqt!(s, λ)  = @timeit to "∇²f" FDDq!(s, λ)
FD2qt!(s, λ)  = @timeit to "∇²f" FD2q!(s, λ)

# Dictionary to hold the results
# Load it if it exists, otherwise creat a new one
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
file_name = joinpath(path_to_package_root, "data", "TimerOutputs")
methods_TimerOutputs_data = Dict()

# make it into short code by listing the methods differently and
# interpolating them using the $ sign
list_timed_methods = [
    ("FD1"        , :qt!,   :Dqt!,  :FDDqt!)
    ("CSD"        , :qt!,   :Dqt!, :CSDDqt!)
    ("DUAL"       , :qt!,   :Dqt!,  :ADDqt!)
    ("FLASH"      , :qt!,   :Dqt!,   :D2qt!)
    ("HYPER"      , :qt!,  :ADqt!,  :AD2qt!)
    ("HYPERSMART" , :qt!, :HSDqt!, :HSD2qt!)
    ("FD2"        , :qt!,  :FDqt!,  :FD2qt!)
]

myruns = ["Compiling run", "Precompiled run"]
myruns2 = ["_compiling_run", ""]

# TimerDict
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

for (i, myrun) in enumerate(myruns)
    for (method_name, q, Dq, D2q) in list_timed_methods
        println("\n\n\n------------------------\n") # print what you are doing
        println(myrun * ": " * method_name * " (with TimerOutputs)")
        println("\n------------------------\n\n\n")

        init.x, init.p = 1x₀, 3p₀ # Reset initial x and p
        J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀ # and J

        # Open the file to print TimerOutputs
        file_name = joinpath(path_to_package_root, "data", "TimerOutputs_" * method_name * myruns2[i] * str_out * ".txt")
        io = open(file_name, "w")

        # Run the timed optimization!
        reset_timer!(to) # Reset timer
        @timeit to "Trust Region" eval(:(res = optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)))

        # Print the TimerOutput
        print_timer(io, to)
        print(io, method_name * "\n" * read(io, String))

        # Close the file
        close(io)

        # Save into Dict
        print_mytimer!(timers, method_name * myruns2[i], to)

    end
end

# Save the results in a Julia data file
using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data", "TimerOutputs_data" * str_out * ".jld2")
@save jld_file timers list_timed_methods



