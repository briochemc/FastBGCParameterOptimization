# Cost functions∇
# Must import the functions to which I add methods
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
    return sqrt(vnorm²(DIN) + vnorm²(POM))
end
function nrm(x::Vector{Dual{U}}) where U
    DIN, POM = unpackx(x)
    return sqrt(vnorm²(DualNumbers.realpart.(DIN)) + vnorm²(DualNumbers.realpart.(POM)))
end
function nrm(x::Vector{Complex{U}}) where U
    DIN, POM = unpackx(x)
    return sqrt(vnorm²(real.(DIN)) + vnorm²(real.(POM)))
end
function nrm(x::Vector{Hyper{U}}) where U
    DIN, POM = unpackx(x)
    return sqrt(vnorm²(HyperDualNumbers.realpart.(DIN)) + vnorm²(HyperDualNumbers.realpart.(POM)))
end

"""
    f(x)

Returns the cost of state `x`.
"""
function f(x) # with respect to x
    DIN, _ = unpackx(x)
    return vnorm²(DIN - DINobs) / vnorm²(DINobs)
end

"""
    ∇ₓf(x)

Returns the gradient of cost of `x` (at `x`).
"""
function ∇ₓf(x)
    DIN, _ = unpackx(x)
    kron([1 0], Dvnorm²(DIN - DINobs) / vnorm²(DINobs))
end
∇ₓf(x, p) = ∇ₓf(x) # for generic form

c_noweight(p::Para) = 0.5 * p2λ(p)' * Matrix(Diagonal(σ²obs.^-1)) * p2λ(p)
c_noweight(p::Para{Complex{Float64}}) = 0.5 * transpose(p2λ(p)) * Matrix(Diagonal(σ²obs.^-1)) * p2λ(p)
Dc_noweight(p::Para) = (∇p2λ(p) .* σ²obs.^-1 .* p2λ(p))'
Dc_noweight(p::Para{Complex{Float64}}) = transpose(∇p2λ(p) .* σ²obs.^-1 .* p2λ(p))

"""
    f(p)

Returns the cost of parameters `p`.
"""
f(p::Para) = p.ω * c_noweight(p)

"""
    ∇ₚf(p)

Returns the gradient of cost of parameters `p` (at `p`).
"""
∇ₚf(p::Para) = p.ω * Dc_noweight(p)
∇ₚf(x, p) = ∇ₚf(p) # for generic form

"""
    f(x, p)

Returns the cost of state `x` plus the cost of parameters `p`.
The costs are added to be used in a Bayesian framework eventually.
(And also because it is simpler.)
"""
function f(x, p::Para) # with respect to both x and p
    return f(x) + f(p)
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
printRMS(cval::Dual) = @printf("RMS = %.2f%% (ε part:%.2g)\n", 100 * sqrt(DualNumbers.realpart(cval) / f(0*x₀)), dualpart(cval))
printRMS(cval::Hyper) = @printf("RMS = %.2f%% (ε₁:%.2g, ε₂:%.2g, ε₁ε₂:%.2g)\n", 100 * sqrt(HyperDualNumbers.realpart(cval) / f(0*x₀)), ε₁part(cval), ε₂part(cval), ε₁ε₂part(cval))
printRMS(cval::Complex) = @printf("RMS = %.2f%% (im part:%.2g)\n", 100 * sqrt(real(cval) / f(0*x₀)), imag(cval))
