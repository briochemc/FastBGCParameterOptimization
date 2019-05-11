
#=============================================
    Constants
=============================================#

"""
    x₀ :: Vector{Float64}

Constant state used to start with.
"""
const nt = 2
const n = nt * nb
const x₀ = ones(n)

"""
    τstop

Constant value for the solver stopping criteria.
Currently set at 1 million years.
"""
const τstop = ustrip(upreferred(1u"Myr"))

#=============================================
    Preallocate buffers for each method
=============================================#

# Preallocate special buffer for F-1 method
println("  Initializing...")
AF1_mem = F1.initialize_mem(x₀, p₀)
println("    AF1 buffer")
F1_mem = F1.initialize_mem(x₀, p₀)
println("    F1 buffer")

# TODO refactor all the methods out of here
list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
for m in list_methods[3:end] # All methods but F1
    mSol = Symbol(string(m) * "Sol")
    @eval $mSol = $m.Solution(copy(x₀))
    println("    $m steady-state solution")
end

str_out = ""

println("Constants and buffers are set up.")
