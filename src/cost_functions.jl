# Cost functions∇
# Must import the functions to which I add methods
import TransportMatrixTools.f̂!, TransportMatrixTools.∇f̂!, TransportMatrixTools.∇²f̂!
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
    f̂!(p; preprint)

Objective `f(sol(p), p)` at `p`.
`F(x, p) = 0` will be solved for a solution `sol` if required.
"""
f̂!(p::Para{Float64}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, init, p, τstop; preprint=preprint)
f̂!(p::Para{Dual{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, εsol, init, p, τstop; preprint=preprint)
HSf̂!(p::Para{Hyper{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, J, hsol, init, p, τstop; preprint=preprint)
f̂!(p::Para{Hyper{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, hsol, init, p, τstop; preprint=preprint)
f̂!(p::Para{Complex{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, imsol, init, p, τstop; preprint=preprint)

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

"""
    f̂!(λ; preprint)

Objective `f(sol(λ), λ)` at `λ`.
`F(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
f̂!(λ::Vector; preprint=" ") = f̂!(λ2p(λ); preprint=preprint)
HSf̂!(λ::Vector; preprint=" ") = HSf̂!(λ2p(λ); preprint=preprint)

"""
    ∇f̂!(λ; preprint)

Evaluates the Gradient of the objective at `λ`.
`F(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
function ∇f̂!(λ::Vector{Float64}; preprint=" ")
    return ∇f̂!(∇ₓf, ∇ₚf, F, ∇ₓF, ∇ₚF, nrm, λ2p, ∇λ2p, J, init, λ, τstop; preprint=preprint)
end
function ∇f̂!(ελ::Vector{Dual{Float64}}; preprint=" ")
    return ∇f̂!(∇ₓf, ∇ₚf, F, ∇ₓF, ∇ₚF, nrm, λ2p, ∇λ2p, εJ, εsol, init, ελ, τstop; preprint=preprint)
end
function ∇f̂!(imλ::Vector{Complex{Float64}}; preprint=" ")
    return ∇f̂!(∇ₓf, ∇ₚf, F, ∇ₓF, ∇ₚF, nrm, λ2p, ∇λ2p, imJ, imsol, init, imλ, τstop; preprint=preprint)
end

"""
    HS∇f̂!(λ; preprint)

Evaluates the HYPERSMART Gradient of the objective at `λ`.
"""
HS∇f̂!(λ::Vector{Float64}; preprint=" ") = ∇f̂!(f, F, ∇ₓF, nrm, λ2p, ∇λ2p, HSbuf, J, hsol, init, λ, τstop; preprint=preprint)

"""
    HS∇²f̂!(λ; preprint)

Evaluates the HYPERSMART Hessian of the objective at `λ`.
"""
HS∇²f̂!(λ::Vector{Float64}; preprint=" ") = ∇²f̂!(f, F, ∇ₓF, nrm, λ2p, ∇λ2p, ∇²λ2p, HSbuf, J, hsol, init, λ, τstop; preprint=preprint)

"""
    ∇²f̂!(λ; preprint)

Evaluates the Hessian of the objective at `λ`.
`F(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
∇²f̂!(λ::Vector{Float64}; preprint=" ") = ∇²f̂!(∇ₓf, ∇ₚf, F, ∇ₓF, ∇ₚF, nrm, λ2p, ∇λ2p, ∇²λ2p, J, init, λ, τstop; preprint=preprint)

"""
    gradient_f̂!(λ; preprint)

Evaluates the gradient of the objective at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
gradient_f̂!(λ::Vector{Float64}) = Calculus.gradient(f̂!, λ)'

"""
    FDf̂!(λ; preprint)

Evaluates the gradient of the objective at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
FDf̂!(λ::Vector{Float64}) = Calculus.jacobian(λ -> [f̂!(λ)], λ, :central)

"""
    FD²f̂!(λ; preprint)

Evaluates the Hessian of the objective at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
FD²f̂!(λ::Vector{Float64}) = Calculus.hessian(f̂!, λ)

"""
    FD∇f̂!(λ; preprint)

Evaluates the Hessian of the objective at `λ` using the `Calculus` package.
The difference with `hessian_q` is that it uses my fast `∇f̂!` and calculates its jacobian.
(Warning: could be slow and inacurrate!)
"""
FD∇f̂!(λ::Vector{Float64}) = Calculus.jacobian(λ -> vec(∇f̂!(λ)), λ, :central)

"""
    CSDf̂!(λ)

Evaluates the gradient of the objective at `λ` using the complex-step method.
(Warning: could be slow!)
"""
CSDf̂!(λ) = ComplexStepGradient(f̂!, λ)'

"""
    CSD∇f̂!(λ)

Evaluates the Hessian of the objective at `λ` using the complex-step method.
Uses the analytical `∇f̂!`.
(Warning: could be slow!)
"""
CSD∇f̂!(λ) = ComplexStepJacobian(∇f̂!, λ)

"""
    ADf̂!(λ)

Evaluates the gradient of the objective at `λ` using DualNumbers.
"""
ADf̂!(λ) = DualNumbersGradient(f̂!, λ)'

"""
    AD∇f̂!(λ)

Evaluates the Hessian of the objective at `λ` using DualNumbers.
Uses the analytical `∇f̂!`.
"""
AD∇f̂!(λ) = DualNumbersJacobian(∇f̂!, λ)

"""
AD²f̂!(λ)

Evaluates the Hessian of the objective at `λ` using HyperDualNumbers.
"""
AD²f̂!(λ) = HyperDualNumbersHessian(f̂!, λ)

"""
    ∇f̂!(storage, λ)

Analytical gradient
"""
function ∇f̂!(storage, λ)
    storage[1:npopt] .= vec(∇f̂!(λ))
end

"""
    ∇²f̂!(storage, λ)

FLASH Hessian
"""
function ∇²f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= ∇²f̂!(λ)
end

"""
    CSDf̂!(storage, λ)

Complex-step method gradient
"""
function CSDf̂!(storage, λ)
    storage[1:npopt] .= vec(CSDf̂!(λ))
end

"""
    FDf̂!(storage, λ)

Finite-difference method gradient
"""
function FDf̂!(storage, λ)
    storage[1:npopt] .= vec(FDf̂!(λ))
end

"""
    ADf̂!(storage, λ)

Dual-step method gradient
"""
function ADf̂!(storage, λ)
    storage[1:npopt] .= vec(ADf̂!(λ))
end

"""
    HS∇f̂!(storage, λ)

HYPERSMART method gradient
"""
function HS∇f̂!(storage, λ)
    storage[1:npopt] .= vec(HS∇f̂!(λ))
end

"""
    CSD∇f̂!(storage, λ)

Complex-step method Hessian
"""
function CSD∇f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= CSD∇f̂!(λ)
end

"""
    FD²f̂!(storage, λ)

Finite-difference method Hessian (2nd order)
"""
function FD²f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= FD²f̂!(λ)
end

"""
    FD∇f̂!(storage, λ)

Finite-difference method Hessian (1st order from analytical gradient)
"""
function FD∇f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= FD∇f̂!(λ)
end

"""
    AD²f̂!(storage, λ)

Hyperdual-step method Hessian
"""
function AD²f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= AD²f̂!(λ)
end

"""
    AD∇f̂!(storage, λ)

Dual-step method Hessian (from the analytical gradient)
"""
function AD∇f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= AD∇f̂!(λ)
end

"""
    HS∇²f̂!(storage, λ)

HYPERSMART method Hessian
"""
function HS∇²f̂!(storage, λ)
    storage[1:npopt, 1:npopt] .= HS∇²f̂!(λ)
end


