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

Objective `f(s(p), p)` at `p`.
`F(x, p) = 0` will be solved for a solution `s(p)` if required.
"""
f̂!(p::Para{Float64}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, init, p, τstop; preprint=preprint)
f̂!(p::Para{Dual{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, εsol, init, p, τstop; preprint=preprint)
F1_f̂!(p::Para{Hyper{Float64}}; preprint=" ") = f̂!(f, F, ∇ₓF, nrm, J, hsol, init, p, τstop; preprint=preprint)
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
F1_f̂!(λ::Vector; preprint=" ") = F1f̂!(λ2p(λ); preprint=preprint)

"""
    ∇f̂!(λ; preprint)

Analytical gradient of the objective at `λ`.
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
    F0_∇²f̂!(λ; preprint)

F0-method Hessian of the objective at `λ`.
"""
F0_∇²f̂!(λ::Vector{Float64}; preprint=" ") = ∇²f̂!(∇ₓf, ∇ₚf, F, ∇ₓF, ∇ₚF, nrm, λ2p, ∇λ2p, ∇²λ2p, J, init, λ, τstop; preprint=preprint)



"""
    F1_∇f̂!(λ; preprint)

F1-method gradient of the objective at `λ`.
"""
F1_∇f̂!(λ::Vector{Float64}; preprint=" ") = ∇f̂!(f, F, ∇ₓF, nrm, λ2p, ∇λ2p, F1buf, J, hsol, init, λ, τstop; preprint=preprint)

"""
    F1_∇²f̂!(λ; preprint)

F1-method Hessian of the objective at `λ`.
"""
F1_∇²f̂!(λ::Vector{Float64}; preprint=" ") = ∇²f̂!(f, F, ∇ₓF, nrm, λ2p, ∇λ2p, ∇²λ2p, F1buf, J, hsol, init, λ, τstop; preprint=preprint)


"""
    OF1_∇f̂!(λ; preprint)

OF1-method gradient of the objective at `λ`.
"""
OF1_∇f̂!(λ::Vector{Float64}; preprint=" ") = ∇f̂!(f, F, ∇ₓf, ∇ₓF, nrm, λ2p, ∇λ2p, ∇sbuf, J, init, λ, τstop; preprint=preprint)

"""
    OF1_∇²f̂!(λ; preprint)

OF1-method Hessian of the objective at `λ`.
"""
OF1_∇²f̂!(λ::Vector{Float64}; preprint=" ") = ∇²f̂!(f, F, ∇ₓf, ∇ₓF, nrm, λ2p, ∇λ2p, ∇²λ2p, ∇sbuf, J, init, λ, τstop; preprint=preprint)

"""
    gradient_f̂!(λ; preprint)

Evaluates the gradient of the objective at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
gradient_f̂!(λ::Vector{Float64}) = Calculus.gradient(f̂!, λ)'

"""
    FD2_∇f̂!(λ; preprint)

FD2-method gradient of the objective at `λ`.
(Warning: could be slow and inacurrate!)
"""
FD2_∇f̂!(λ::Vector{Float64}) = Calculus.jacobian(λ -> [f̂!(λ)], λ, :central)

"""
    FD2_∇²f̂!(λ; preprint)

FD2-method Hessian of the objective at `λ`.
(Warning: could be slow and inacurrate!)
"""
FD2_∇²f̂!(λ::Vector{Float64}) = Calculus.hessian(f̂!, λ)

"""
    FD1_∇²f̂!(λ; preprint)


FD1-method Hessian of the objective at `λ`.
"""
FD1_∇²f̂!(λ::Vector{Float64}) = Calculus.jacobian(λ -> vec(∇f̂!(λ)), λ, :central)

"""
    CSD_∇f̂!(λ)

CSD method used to compute the gradient of the objective at `λ`.
(Not part of the CSD method in the paper, which uses the analytical gradient.)
"""
CSD_∇f̂!(λ) = ComplexStepGradient(f̂!, λ)'

"""
    CSD_∇²f̂!(λ)

CSD-method Hessian of the objective at `λ`.
"""
CSD_∇²f̂!(λ) = ComplexStepJacobian(∇f̂!, λ)

"""
    DUAL_∇²f̂!(λ)

DUAL-method Hessian of the objective at `λ`.
"""
DUAL_∇²f̂!(λ) = DualNumbersJacobian(∇f̂!, λ)

"""
    HYPER_∇f̂!(λ)

HYPER-method gradient of the objective at `λ`.
(Note HYPER is naive and uses dual numbers here.)
"""
HYPER_∇f̂!(λ) = DualNumbersGradient(f̂!, λ)'

"""
    HYPER_∇²f̂!(λ)

HYPER-method Hessian of the objective at `λ`.
"""
HYPER_∇²f̂!(λ) = HyperDualNumbersHessian(f̂!, λ)

#=
    Objective and derivatives with storage
=#

# TODO shorten this all in an expression loop

Hessian_methods = [
    :OF1_∇²f̂!,
    :F1_∇²f̂!,
    :F0_∇²f̂!,
    :DUAL_∇²f̂!,
    :CSD_∇²f̂!,
    :HYPER_∇²f̂!,
    :FD1_∇²f̂!,
    :FD2_∇²f̂!
]

gradient_methods = [
    :∇f̂!,
    :OF1_∇f̂!,
    :F1_∇f̂!,
    :CSD_∇f̂!,
    :HYPER_∇f̂!,
    :FD2_∇f̂!
]

for m in Hessian_methods
    @eval $m(s, λ) = s[1:npopt, 1:npopt] .= $m(λ)
end

for m in gradient_methods
    @eval $m(s, λ) = s[1:npopt] .= vec($m(λ))
end

#"""
#    ∇f̂!(storage, λ)
#
#Analytical gradient
#"""
#function ∇f̂!(storage, λ)
#    storage[1:npopt] .= vec(∇f̂!(λ))
#end
#
#"""
#    F0_∇²f̂!(storage, λ)
#
#F0 Hessian
#"""
#function F0_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= F0_∇²f̂!(λ)
#end
#
#"""
#    CSD_∇²f̂!(storage, λ)
#
#CSD Hessian
#"""
#function CSD_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= CSD_∇²f̂!(λ)
#end
#
#"""
#    DUAL_∇²f̂!(storage, λ)
#
#DUAL Hessian
#"""
#function DUAL_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= DUAL_∇²f̂!(λ)
#end
#
#"""
#    FD1_∇²f̂!(storage, λ)
#
#FD1 Hessian
#"""
#function FD1_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= FD1_∇²f̂!(λ)
#end
#
#"""
#    F1_∇f̂!(storage, λ)
#
#F1 gradient
#"""
#function F1_∇f̂!(storage, λ)
#    storage[1:npopt] .= vec(F1_∇f̂!(λ))
#end
#
#"""
#    F1_∇²f̂!(storage, λ)
#
#F1 Hessian
#"""
#function F1_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= F1_∇²f̂!(λ)
#end
#
#"""
#    FD2_∇f̂!(storage, λ)
#
#FD2 gradient
#"""
#function FD2_∇f̂!(storage, λ)
#    storage[1:npopt] .= vec(FD2_∇f̂!(λ))
#end
#
#"""
#    FD2_∇²f̂!(storage, λ)
#
#FD2 Hessian
#"""
#function FD2_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= FD2_∇²f̂!(λ)
#end
#
#"""
#    HYPER_∇f̂!(storage, λ)
#
#HYPER gradient
#"""
#function HYPER_∇f̂!(storage, λ)
#    storage[1:npopt] .= vec(HYPER_∇f̂!(λ))
#end
#
#"""
#    HYPER_∇²f̂!(storage, λ)
#
#HYPER Hessian
#"""
#function HYPER_∇²f̂!(storage, λ)
#    storage[1:npopt, 1:npopt] .= HYPER_∇²f̂!(λ)
#end
#
