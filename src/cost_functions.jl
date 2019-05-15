
#=======
using AIBECS
=======#
# Hyper parameters
ωs = [1.0, 0.0, 0.0]
ωp = 1e-6
# PO₄ mean and variance of observations fom WOA18
using WorldOceanAtlasTools
WOA = WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WOA.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
const μx = (μDIPobs, missing, missing)
const σ²x = (σ²DIPobs, missing, missing)
# generate mismatch functions
f = generate_objective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)
∇ₓf = generate_∇ₓobjective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)
∇ₚf = generate_∇ₚobjective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)



# Cost functions
# TODO code inside AIBECS?
"""
    nrm(x)

Gives the "tracer" norm:
`nrm(x)` is the square root of the sum the squares of the volume-weighted norms of (the real parts of) each tracer.
This choice is arbitrary but is simple.
However, it has a weird unit of mol m^(-3/2).
This is OK because we normalize everything later.
(Particularly important for the cost function and the solver tolerances.)
Note: `nrm` **should not** be used with the complex step method or dual numbers.
"""
function nrm(x)
    DIP, DOP, POP = unpackx(x)
    return sqrt(vnorm²(ℜ(DIP)) + vnorm²(ℜ(DOP)) + vnorm²(ℜ(POP)))
end
ℜ(x::Real) = (x)
ℜ(x::Complex) = real(x)
ℜ(x::Dual) = DualNumbers.realpart(x)
ℜ(x::Hyper) = HyperDualNumbers.realpart(x)
ℜ(x::Vector) = ℜ.(x)

vnorm²(x) = AIBECS.weighted_norm²(x, v)
"""
    print_cost(cval; preprint)

Prints the cost as a root mean square (RMS) error in percent.
(Will also print the imaginary or dual part if any.)
"""
function print_cost(cval; preprint=" ")
    if preprint ≠ ""
        print(preprint)
        printRMS(cval)
    end
    return nothing
end
printRMS(cval) = @printf("RMS = %.2f%%\n", 100 * sqrt(cval / f(0*x₀)))
printRMS(cval::Dual) = @printf("RMS = %.2f%% (ε part:%.2g)\n", 100 * sqrt(ℜ(cval) / f(0*x₀)), 𝔇(cval))
printRMS(cval::Hyper) = @printf("RMS = %.2f%% (ε₁:%.2g, ε₂:%.2g, ε₁ε₂:%.2g)\n", 100 * sqrt(ℌ(cval) / f(0*x₀)), ℌ₁(cval), ℌ₂(cval), ℌ(cval))
printRMS(cval::Complex) = @printf("RMS = %.2f%% (im part:%.2g)\n", 100 * sqrt(ℜ(cval) / f(0*x₀)), ℑ(cval))

ℑ(x::Complex) = imag(x)
𝔇(x::Dual) = DualNumbers.dualpart(x)
ℌ(x::Dual) = HyperDualNumbers.ε₁ε₂part(x)
ℌ₁(x::Dual) = HyperDualNumbers.ε₁part(x)
ℌ₂(x::Dual) = HyperDualNumbers.ε₂part(x)
