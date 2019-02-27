"""
    p₀

The (constant) default values of non-optimized parameters.
"""
const p₀ = Para(ω = 0.0)
str_out = "_w0"

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
    τstop

Constant value for the solver stopping criteria.
Currently set at 1 million years.
"""
const τstop = 1e6 * 365e6 * spd

# Preallocate real, dual, complex, and hyperdual states (and Jacobians)
init, J, εsol, εJ, imsol, imJ, hsol = preallocateNewTypes(Para, fJac, x₀, p₀)


