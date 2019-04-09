
#=============================================
    Constants
=============================================#

"""
    x₀ :: Vector{Float64}

Constant state used to start with.
"""
const x₀ = [DINobs; DINobs / 10]

"""
    n :: Int

Constant state used to start with.
"""
const n = length(x₀)

"""
    τstop

Constant value for the solver stopping criteria.
Currently set at 1 million years.
"""
const τstop = 1e6 * 365e6 * spd

#=============================================
    Preallocate buffers for each method
=============================================#

# Preallocate special buffer for F-1 method
println("  Initializing...")
F1buf = F1.initialize_buffer(x₀, p₀)
println("    F1 buffer")

# TODO refactor all the methods out of here
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
for m in list_methods[2:end] # All methods but F1
    mSol = Symbol(string(m) * "Sol")
    @eval $mSol = $m.Solution(copy(x₀))
    println("    $m steady-state solution")
end

println("Constants and buffers are set up.")
