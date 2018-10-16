# Cost functions
import TransportMatrixTools.q!, TransportMatrixTools.Dq!, TransportMatrixTools.D2q!

function nrm(x)
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(DSi) + vnorm²(PSi))
end# This could change if I want with some other norm
function nrm(x::Vector{Dual{U}}) where U
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(realpart.(DSi)) + vnorm²(realpart.(PSi)))
end# This could change if I want with some other norm

# cost function
function c(x) # with respect to x
    DSi, _ = unpackx(x)
    return vnorm²(DSi - DSiobs) / vnorm²(DSiobs)
end
function Dc(x)
    DSi, _ = unpackx(x)
    kron([1 0], Dvnorm²(DSi - DSiobs) / vnorm²(DSiobs))
end

c_noweight(p::Para) = 0.5 * p2λ(p)' * p2λ(p)
Dc_noweight(p::Para) = Dp2λ(p)' .* p2λ(p)'
const x₀ = [DSiobs; DSiobs / 10] * 1.1
const ω = 1e-2 * c(x₀) # To be determined!
c(p::Para) = ω * c_noweight(p)
Dc(p::Para) = ω * Dc_noweight(p)

function c(x, p::Para) # with respect to both x and p
    return c(x) + c(p)
end

function q!(init::RealSolution, p::Para{Float64}, c, f, fJac, nrm, τstop, verbose::Bool)
    if verbose
        show(IOContext(stdout, :compact => true), p)
        return q!(init, p, c, f, nrm, τstop, "    ")
    else
        return q!(init, p, c, f, nrm, τstop, "")
    end
end
function q!(εsol::DualSolution, init::RealSolution, εp::Para{Dual{Float64}}, c, f, fJac, nrm, τstop, verbose::Bool)
    if verbose
        show(IOContext(stdout, :compact => true), εp)
        return q!(εsol, init, εp, c, f, nrm, τstop, "    ")
    else
        return q!(εsol, init, εp, c, f, nrm, τstop, "")
    end
end

# Preallocate initial state and Jacobian, and τstop for wrapping qprint
const λ₀ = p2λ(p₀)
# SaJ = StateAndJacobian(x₀, factorize(fJac(x₀, p₀)), fJac(x₀, p₀), p₀) # the Jacobian factors
init = RealSolution(x₀, p₀)
MyJ = MyRealJacobianFactors(factorize(fJac(x₀, p₀)), p₀)
# Also preallocate the Dual containers
εp₀ = Para(convert(Vector{Dual{Float64}}, vec(p₀)))
εx₀ = convert(Vector{Dual{Float64}}, x₀)
εsol = DualSolution(εx₀, εp₀)
MyεJ = MyDualJacobianFactors(factorize(fJac(εx₀, εp₀)), εp₀)
const τstop = 1e6 * 365e6 * spd
q!(p::Para{Float64}) = q!(init, p, c, f, fJac, nrm, τstop, false)
q!(p::Para{Dual{Float64}}) = q!(εsol, init, p, c, f, fJac, nrm, τstop, false)

# Need to define the function with a storage argument first
function q!(init::RealSolution, λ::Vector{Float64}, c, f, fJac, nrm, τstop, λ2p, verbose::Bool)
    if verbose
        show(IOContext(stdout, :compact => true), λ2p(λ))
        qval = q!(init, λ, c, f, fJac, nrm, τstop, λ2p, "    ")
        print_cost(qval, "    ")
        return qval
    else
        qval = q!(init, λ, c, f, fJac, nrm, τstop, λ2p, "")
        return qval
    end
end
function q!(εsol::DualSolution, init::RealSolution, ελ::Vector{Dual{Float64}}, c, f, fJac, nrm, τstop, λ2p, verbose::Bool)
    if verbose
        show(IOContext(stdout, :compact => true), λ2p(ελ))
        qval = q!(εsol, init, ελ, c, f, fJac, nrm, τstop, λ2p, "    ")
        print_cost(qval, "    ")
        return qval
    else
        qval = q!(εsol, init, ελ, c, f, fJac, nrm, τstop, λ2p, "")
        return qval
    end
end

function print_cost(cval, preprint = "")
    if preprint ≠ ""
        print(preprint)
        @printf("RMS = %.2f%%\n", 100 * sqrt(cval / c(0*x₀)))
    end
    return nothing
end
function print_cost(cval::Dual{Float64}, preprint = "")
    if preprint ≠ ""
        print(preprint)
        @printf("RMS = %.2f%% (+ ε part:%.2g)\n", 100 * sqrt(realpart(cval) / c(0*x₀)), dualpart(cval))
    end
    return nothing
end

q!(λ::Vector{Float64}) = q!(init, λ, c, f, fJac, nrm, τstop, λ2p, false)
q!(ελ::Vector{Dual{Float64}}) = q!(εsol, init, ελ, c, f, fJac, nrm, τstop, λ2p, true)
# Qwrap(λ) = Q!(c, λ, SaJ, f, Dxf, vnorm, τstop)
# slowQwrap(λ) = slowQ(c, λ, nwet, f, fJac, nrm, τstop)
# vnorm(f(ones(length(x₀)),p₀))
# Q₀ = Qwrap(λ₀)
Dq!(λ::Vector{Float64}) = Dq!(MyJ, init, λ, Dc, f, fJac, Dpf, nrm, τstop, λ2p, Dλ2p, "")
Dq!(ελ::Vector{Dual{Float64}}) = Dq!(MyεJ, εsol, init, ελ, Dc, f, fJac, Dpf, nrm, τstop, λ2p, Dλ2p, "")
# slowDQwrap(λ) = slowDQ(Dc, λ, nwet, f, fJac, Dpf, nrm, τstop)
D2q!(λ::Vector{Float64}) = D2q!(MyJ, init, λ, Dc, f, fJac, Dpf, nrm, τstop, λ2p, Dλ2p, D2λ2p, "")

slowDq!(λ::Vector{Float64}) = Calculus.gradient(q!, λ)'
slowD2q!(λ::Vector{Float64}) = Calculus.hessian(q!, λ)


function Dq!(storage, λ)
    storage[1:npopt] .= vec(Dq!(λ))
end
function D2q!(storage, λ)
    storage[1:npopt, 1:npopt] .= D2q!(λ)
end

function slowDq!(storage, λ)
    storage[1:npopt] .= vec(slowDq!(λ))
end
function slowD2q!(storage, λ)
    storage[1:npopt, 1:npopt] .= slowD2q!(λ)
end
