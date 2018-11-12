# Cost functions

# Must import the functions to which I add methods
import TransportMatrixTools.q!, TransportMatrixTools.Dq!, TransportMatrixTools.D2q!

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
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(DSi) + vnorm²(PSi))
end
function nrm(x::Vector{Dual{U}}) where U
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(realpart.(DSi)) + vnorm²(realpart.(PSi)))
end
function nrm(x::Vector{Complex{U}}) where U
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(real.(DSi)) + vnorm²(real.(PSi)))
end

"""
    c(x)

Returns the cost of state `x`.
"""
function c(x) # with respect to x
    DSi, _ = unpackx(x)
    return vnorm²(DSi - DSiobs) / vnorm²(DSiobs)
end

"""
    Dc(x)

Returns the gradient of cost of `x` (at `x`).
"""
function Dc(x)
    DSi, _ = unpackx(x)
    kron([1 0], Dvnorm²(DSi - DSiobs) / vnorm²(DSiobs))
end

c_noweight(p::Para) = 0.5 * p2λ(p)' * p2λ(p)
c_noweight(p::Para{Complex{Float64}}) = 0.5 * transpose(p2λ(p)) * p2λ(p)
Dc_noweight(p::Para) = Dp2λ(p)' .* p2λ(p)'
Dc_noweight(p::Para{Complex{Float64}}) = transpose(Dp2λ(p) .* p2λ(p))

#c(p::Para) = 0.5 * 
#δ(p::Para) = p2λ(p)



"""
    x₀ :: Vector{Float64}

Constant state used to start with.
"""
const x₀ = [DSiobs; DSiobs / 10] * 1.1

"""
    ω :: Float64

Constant parameter for the relative weight of the cost of `p` relative to the cost of `x`.
Can be used for Bayesian priors I guess.
Right now used to impose a parabolic shape to the overall cost function (avoids compensations).
"""
const ω = 1e-2 * c(x₀) # To be determined!

"""
    c(p)

Returns the cost of parameters `p`.
"""
c(p::Para) = ω * c_noweight(p)

"""
    Dc(p)

Returns the gradient of cost of parameters `p` (at `p`).
"""
Dc(p::Para) = ω * Dc_noweight(p)

"""
    c(x, p)

Returns the cost of state `x` plus the cost of parameters `p`.
The costs are added to be used in a Bayesian framework eventually.
(And also because it is simpler.)
"""
function c(x, p::Para) # with respect to both x and p
    return c(x) + c(p)
end

# Preallocate real, dual, complex, and hyperdual states (and Jacobians)
init, J, εsol, εJ, imsol, imJ, hsol = preallocateNewTypes(Para, fJac, x₀, p₀)

"""
    τstop

Constant value for the solver stopping criteria.
Currently set at 1 million years.
"""
const τstop = 1e6 * 365e6 * spd

"""
    q!(p; preprint)

Full cost `c(sol(p), p)` at `p`.
`f(x, p) = 0` will be solved for a solution `sol` if required.
"""
q!(p::Para{Float64}; preprint="") = q!(c, f, fJac, nrm, init, p, τstop; preprint=preprint)
q!(p::Para{Dual{Float64}}; preprint="") = q!(c, f, fJac, nrm, εsol, init, p, τstop; preprint=preprint)
q!(p::Para{Hyper{Float64}}; preprint="") = q!(c, f, fJac, nrm, hsol, init, p, τstop; preprint=preprint)
q!(p::Para{Complex{Float64}}; preprint="") = q!(c, f, fJac, nrm, imsol, init, p, τstop; preprint=preprint)

"""
    print_cost(cval; preprint)

Prints the cost as a root mean square (RMS) error in percent.
(Will also print the imaginary or dual part if any.)
"""
function print_cost(cval; preprint = "")
    if preprint ≠ ""
        print(preprint)
        printRMS(cval)
    end
    return nothing
end
printRMS(cval) = @printf("RMS = %.2f%%\n", 100 * sqrt(cval / c(0*x₀)))
printRMS(cval::Dual) = @printf("RMS = %.2f%% (ε part:%.2g)\n", 100 * sqrt(realpart(cval) / c(0*x₀)), dualpart(cval))
printRMS(cval::Hyper) = @printf("RMS = %.2f%% (ε₁:%.2g, ε₂:%.2g, ε₁ε₂:%.2g)\n", 100 * sqrt(real(cval) / c(0*x₀)), eps1(cval), eps2(cval), eps1eps2(cval))
printRMS(cval::Complex) = @printf("RMS = %.2f%% (im part:%.2g)\n", 100 * sqrt(real(cval) / c(0*x₀)), imag(cval))

"""
    q!(λ; preprint)

Full cost `c(sol(λ), λ)` at `λ`.
`f(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
q!(λ::Vector; preprint="") = q!(λ2p(λ); preprint=preprint)

"""
    Dq!(λ; preprint)

Evaluates the Gradient of the full cost at `λ`.
`f(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
function Dq!(λ::Vector{Float64}; preprint="")
    return Dq!(Dc, f, fJac, Dpf, nrm, λ2p, Dλ2p, J, init, λ, τstop; preprint=preprint)
end
function Dq!(ελ::Vector{Dual{Float64}}; preprint="")
    return Dq!(Dc, f, fJac, Dpf, nrm, λ2p, Dλ2p, εJ, εsol, init, ελ, τstop; preprint=preprint)
end
function Dq!(imλ::Vector{Complex{Float64}}; preprint="")
    return Dq!(Dc, f, fJac, Dpf, nrm, λ2p, Dλ2p, imJ, imsol, init, imλ, τstop; preprint=preprint)
end

"""
    D2q!(λ; preprint)

Evaluates the Hessian of the full cost at `λ`.
`f(x, p(λ)) = 0` will be solved for a solution `sol` if required.
"""
D2q!(λ::Vector{Float64}; preprint="") = D2q!(Dc, f, fJac, Dpf, nrm, λ2p, Dλ2p, D2λ2p, J, init, λ, τstop; preprint=preprint)

"""
    gradient_q!(λ; preprint)

Evaluates the gradient of the full cost at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
gradient_q!(λ::Vector{Float64}) = Calculus.gradient(q!, λ)'

"""
    FDq!(λ; preprint)

Evaluates the gradient of the full cost at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
FDq!(λ::Vector{Float64}) = Calculus.jacobian(λ -> [q!(λ)], λ, :central)

"""
    FD2q!(λ; preprint)

Evaluates the Hessian of the full cost at `λ` using the `Calculus` package.
(Warning: could be slow and inacurrate!)
"""
FD2q!(λ::Vector{Float64}) = Calculus.hessian(q!, λ)

"""
    FDDq!(λ; preprint)

Evaluates the Hessian of the full cost at `λ` using the `Calculus` package.
The difference with `hessian_q` is that it uses my fast `Dq!` and calculates its jacobian.
(Warning: could be slow and inacurrate!)
"""
FDDq!(λ::Vector{Float64}) = Calculus.jacobian(λ -> vec(Dq!(λ)), λ, :central)

"""
    ComplexStepGradient(q, λ::Vector{U})

Returns the gradient using the complex step method.
Only good for small sizes.
`q` is an application from Rⁿ to R in this case.
"""
function ComplexStepGradient(q, λ::Vector{U}) where U # q : Rⁿ -> R
    n = length(λ)
    out = zeros(U, n)
    h = 1e-50
    for i in 1:length(λ)
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[i] = imag.(q(imλ)) / h
    end
    return out
end

"""
    CSDq!(λ)

Evaluates the gradient of the full cost at `λ` using the complex-step method.
(Warning: could be slow!)
"""
CSDq!(λ) = ComplexStepGradient(q!, λ)'

"""
    ComplexStepJacobian(Dq, λ::Vector{U})

Returns the Jacobian using the complex step method.
Only good for small sizes.
(Do not use for the state model `J` if using OCIM!)
`Dq` is an application from Rⁿ to Rⁿ in this case.
"""
function ComplexStepJacobian(Dq, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    h = 1e-50
    for i in 1:n
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[:, i] .= imag.(vec(Dq(imλ))) / h
    end
    return out
end
"""
    CSDDq!(λ)

Evaluates the Hessian of the full cost at `λ` using the complex-step method.
Uses the analytical `Dq!`.
(Warning: could be slow!)
"""
CSDDq!(λ) = ComplexStepJacobian(Dq!, λ)

function Dq!(storage, λ)
    storage[1:npopt] .= vec(Dq!(λ))
end
function D2q!(storage, λ)
    storage[1:npopt, 1:npopt] .= D2q!(λ)
end

function CSDq!(storage, λ)
    storage[1:npopt] .= vec(CSDq!(λ))
end
function FDq!(storage, λ)
    storage[1:npopt] .= vec(FDq!(λ))
end

function CSDDq!(storage, λ)
    storage[1:npopt, 1:npopt] .= CSDDq!(λ)
end
function FD2q!(storage, λ)
    storage[1:npopt, 1:npopt] .= FD2q!(λ)
end
function FDDq!(storage, λ)
    storage[1:npopt, 1:npopt] .= FDDq!(λ)
end
