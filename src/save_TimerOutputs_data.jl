# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, extended_trace = true, x_tol = 1e-3)
using TimerOutputs
const to = TimerOutput()

# timed version for all functions called in Optim
qt!(λ)        = @timeit to "q"                  q!(λ)
Dqt!(s, λ)    = @timeit to "Dq"                 Dq!(s, λ)
D2qt!(s, λ)   = @timeit to "D2q (mine)"         D2q!(s, λ)
ADDqt!(s, λ)  = @timeit to "D2q (dual)"         ADDq!(s, λ)
CSDDqt!(s, λ) = @timeit to "D2q (complex)"      CSDDq!(s, λ)
FDDqt!(s, λ)  = @timeit to "D2q (finite-diff)"  FDDq!(s, λ)

# Dictionary to hold the results
# Load it if it exists, otherwise creat a new one
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
file_name = joinpath(path_to_package_root, "data", "TimerOutputs")
methods_TimerOutputs_data = Dict()

# make it into short code by listing the methods differently and
# interpolating them using the $ sign
list_timed_methods = [
    (  "D2q", :qt!, :Dqt!,   :D2qt!)
    ( "ADDq", :qt!, :Dqt!,  :ADDqt!)
    ("CSDDq", :qt!, :Dqt!, :CSDDqt!)
    ( "FDDq", :qt!, :Dqt!,  :FDDqt!)
]

myruns = ["Run 1", "Run 2"]
myruns2 = ["_Run1", "_Run2"]

for (i, myrun) in enumerate(myruns)
    for (method_name, q, Dq, D2q) in list_timed_methods
        println("\n\n\n------------------------\n") # print what you are doing
        println(myrun * ": " * method_name * " (with TimerOutputs)")
        println("\n------------------------\n\n\n")

        init.x, init.p = 1x₀, 3p₀ # Reset initial x and p
        J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀ # and J

        # Open the file to print TimerOutputs
        file_name = joinpath(path_to_package_root, "data", "TimerOutputs_" * method_name * myruns2[i] * ".txt")
        io = open(file_name, "w")

        # Run the timed optimization!
        reset_timer!(to) # Reset timer
        @timeit to "Trust Region" eval(:(res = optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)))

        # Print the TimerOutput
        print_timer(io, to)

        # Close the file
        close(io)

    end
end


