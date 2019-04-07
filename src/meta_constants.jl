"""
    p₀

The (constant) default values of non-optimized parameters.
"""
const p₀ = Para(ω = 1e-4)
str_out = "_default"

"""
    μobs

The (constant) mean of the log of the observed parameters (the μ of the lognormal prior).
"""
const μobs = [meanlogx.(metaflatten(p₀, prior))...]

"""
    σ²obs

The (constant) variance of the log of the observed parameters (the σ² of the lognormal prior).
"""
const σ²obs = [varlogx.(metaflatten(p₀, prior))...]

"""
    λ₀

Default `λ` corresponding to `p₀` (constant).
"""
const λ₀ = p2λ(p₀)

"""
    x₀ :: Vector{Float64}

Constant state used to start with.
"""
const x₀ = [DINobs; DINobs / 10] * 1.1

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

# Preallocate real, dual, complex, and hyperdual states (and Jacobians)
init, J, εsol, εJ, imsol, imJ, hsol, F1buf, ∇sbuf = preallocateNewTypes(Para, ∇ₓF, x₀, p₀)

# Preallocate special buffer for F-1 method
println("Initializing FormulaOne Buffer...")
newF1buf = initialize_buffer(f, F, ∇ₓf, ∇ₓF, x₀, p₀, CTKAlg(), nrm=nrm, preprint=" ")
newF1buf.s = 1x₀ # (`1x₀` is equivalent to `copy(x₀)` I think)
# TODO maybe change from initialize_buffer to simply calling the Buffer constructor
# so that this method does not start with an advantage?
# TODO refactor all the methods out of here

println("Constants are set up.")
