# Cost functions
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
    DIN, POM = unpackx(x)
    return sqrt(vnorm²(ℜ(DIN)) + vnorm²(ℜ(POM)))
end
ℜ(x::Real) = (x)
ℜ(x::Complex) = real(x)
ℜ(x::Dual) = DualNumbers.realpart(x)
ℜ(x::Hyper) = HyperDualNumbers.realpart(x)
ℜ(x::Vector) = ℜ.(x)

"""
    fₓ(x)

Returns the cost of state `x`.
"""
function fₓ(x) # with respect to x
    DIN, _ = unpackx(x)
    return vnorm²(DIN - DINobs) / vnorm²(DINobs)
end

"""
    ∇ₓf(x)

Returns the gradient of cost of `x` (at `x`).
"""
function ∇ₓf(x, p)
    DIN, _ = unpackx(x)
    kron([1 0], Dvnorm²(DIN - DINobs) / vnorm²(DINobs))
end

fₚ_noweight(p) = 0.5 * transpose(p2λ(p)) * Matrix(Diagonal(σ²obs.^-1)) * p2λ(p)
Dfₚ_noweight(p) = transpose(∇p2λ(p) .* σ²obs.^-1 .* p2λ(p))

"""
    f(p)

Returns the cost of parameters `p`.
"""
fₚ(p) = p.ω * fₚ_noweight(p)

"""
    ∇ₚf(p)

Returns the gradient of cost of parameters `p` (at `p`).
"""
∇ₚf(x, p) = p.ω * Dfₚ_noweight(p) # for generic form

"""
    f(x, p)

Returns the cost of state `x` plus the cost of parameters `p`.
The costs are added to be used in a Bayesian framework eventually.
(And also because it is simpler.)
"""
function f(x, p) # with respect to both x and p
    return fₓ(x) + fₚ(p)
end


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
