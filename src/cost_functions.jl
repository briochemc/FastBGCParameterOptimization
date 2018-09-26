# Cost functions
import TransportMatrixtools.q!, TransportMatrixtools.Dq!, TransportMatrixtools.D2q!

function nrm(x)
    DSi, PSi = unpackx(x)
    return sqrt(vnorm²(DSi) + vnorm²(PSi))
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

function q!(c, p::Para, SaJ, f, fJac, nrm, τstop, verbose::Bool)
    if verbose
        paraprint(p, "  ")
        return q!(c, p, SaJ, f, fJac, nrm, τstop, "  ")
    else
        return q!(c, p, SaJ, f, fJac, nrm, τstop, "")
    end
end

# Preallocate State and Jacobian, and τstop for wrapping qprint
const λ₀ = p2λ(p₀)
SaJ = StateAndJacobian(x₀, factorize(fJac(x₀, p₀)), p₀) # the Jacobian factors
const τstop = τg * 1e6
q!(p::Para) = q!(c, p, SaJ, f, fJac, nrm, τstop, false)

# Need to define the function with a storage argument first
function q!(c, λ::Vector, SaJ, f, fJac, nrm, τstop, λ2p, verbose::Bool)
    if verbose
        paraprint(λ2p(λ), "    ")
        qval = q!(c, λ, SaJ, f, fJac, nrm, τstop, λ2p, "    ")
        print_cost(qval, "    ")
        return qval
    else
        qval = q!(c, λ, SaJ, f, fJac, nrm, τstop, λ2p, "")
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

q!(λ::Vector) = q!(c, λ, SaJ, f, fJac, nrm, τstop, λ2p, false)
# Qwrap(λ) = Q!(c, λ, SaJ, f, Dxf, vnorm, τstop)
# slowQwrap(λ) = slowQ(c, λ, nwet, f, fJac, nrm, τstop)
# vnorm(f(ones(length(x₀)),p₀))
# Q₀ = Qwrap(λ₀)
Dq!(λ) = Dq!(Dc, λ, SaJ, f, fJac, Dpf, nrm, τstop, λ2p, Dλ2p, "")
# slowDQwrap(λ) = slowDQ(Dc, λ, nwet, f, fJac, Dpf, nrm, τstop)
D2q!(λ) = D2q!(Dc, λ, SaJ, f, fJac, Dpf, nrm, τstop, λ2p, Dλ2p, D2λ2p, "")


function Dq!(storage, λ)
    storage[1:npopt] .= vec(Dq!(λ))
end
function D2q!(storage, λ)
    storage[1:npopt, 1:npopt] .= D2q!(λ)
end


