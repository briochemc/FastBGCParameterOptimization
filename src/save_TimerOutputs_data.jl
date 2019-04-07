# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = false, show_trace = false, extended_trace = false)
using TimerOutputs
const to = TimerOutput()

# timed version for all functions called in Optim
f̂t!(λ)           = @timeit to "f"   f̂!(λ)
f̂t(λ)            = @timeit to "f"   f̂(λ)
∇f̂t!(s, λ)       = @timeit to "∇f"  ∇f̂!(s, λ)       #┐
∇f̂t(s, λ)        = @timeit to "∇f"  ∇f̂(s, λ)        #│
OF1_∇f̂t!(s, λ)   = @timeit to "∇f"  OF1_∇f̂!(s, λ)   #│
F1_∇f̂t!(s, λ)    = @timeit to "∇f"  F1_∇f̂!(s, λ)    #│
HYPER_∇f̂t!(s, λ) = @timeit to "∇f"  HYPER_∇f̂!(s, λ) #├─ gradients
FD2_∇f̂t!(s, λ)   = @timeit to "∇f"  FD2_∇f̂!(s, λ)   #┘
OF1_∇²f̂t!(s, λ)  = @timeit to "∇²f" OF1_∇²f̂!(s, λ)   #┐
∇²f̂t(s, λ)       = @timeit to "∇²f" ∇²f̂(s, λ)        #│
F0_∇²f̂t!(s, λ)   = @timeit to "∇²f" F0_∇²f̂!(s, λ)    #│
F1_∇²f̂t!(s, λ)   = @timeit to "∇²f" F1_∇²f̂!(s, λ)    #├─ Hessians
DUAL_∇²f̂t!(s, λ) = @timeit to "∇²f" DUAL_∇²f̂!(s, λ)  #│
CSD_∇²f̂t!(s, λ)  = @timeit to "∇²f" CSD_∇²f̂!(s, λ)   #│
FD1_∇²f̂t!(s, λ)  = @timeit to "∇²f" FD1_∇²f̂!(s, λ)   #│
HYPER_∇²f̂t!(s, λ)= @timeit to "∇²f" HYPER_∇²f̂!(s, λ) #│
FD2_∇²f̂t!(s, λ)  = @timeit to "∇²f" FD2_∇²f̂!(s, λ)   #┘

# Dictionary to hold the results
# Load it if it exists, otherwise creat a new one
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
file_name = joinpath(path_to_package_root, "data", "TimerOutputs")
methods_TimerOutputs_data = Dict()

# make it into short code by listing the methods differently and
# interpolating them using the $ sign
list_timed_methods = [
    ("newF1" ,  :f̂t,       :∇f̂t!,        :∇²f̂t)
    ("OF1"   , :f̂t!,   :OF1_∇f̂t!,   :OF1_∇²f̂t!)
    ("F0"    , :f̂t!,       :∇f̂t!,    :F0_∇²f̂t!)
    ("F1"    , :f̂t!,       :∇f̂t!,    :F1_∇²f̂t!)
    ("DUAL"  , :f̂t!,       :∇f̂t!,  :DUAL_∇²f̂t!)
    ("CSD"   , :f̂t!,       :∇f̂t!,   :CSD_∇²f̂t!)
    ("FD1"   , :f̂t!,       :∇f̂t!,   :FD1_∇²f̂t!)
    ("HYPER" , :f̂t!, :HYPER_∇f̂t!, :HYPER_∇²f̂t!)
    ("FD2"   , :f̂t!,   :FD2_∇f̂t!,   :FD2_∇²f̂t!)
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
    for (method_name, objective, gradient, hessian) in list_timed_methods
        println("\n\n\n------------------------\n") # print what you are doing
        println(myrun * ": " * method_name * " (with TimerOutputs)")
        println("\n------------------------\n\n\n")

        init.x, init.p = 1x₀, 3p₀ # Reset initial x and p
        J.fac, J.p = factorize(∇ₓF(x₀, 3p₀)), 3p₀ # and J
        newF1buf.s, newF1buf.p = 1x₀, 3p₀

        # Open the file to print TimerOutputs
        file_name = joinpath(path_to_package_root, "data", "TimerOutputs_" * method_name * myruns2[i] * str_out * ".txt")
        io = open(file_name, "w")

        # Run the timed optimization!
        reset_timer!(to) # Reset timer
        @timeit to "Trust Region" eval(:(res = optimize($objective, $gradient, $hessian, $λ₀, NewtonTrustRegion(), $opt)))

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



