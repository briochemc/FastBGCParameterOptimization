# Geological restoring
function geores(DIN, p::Para)
    τg = p.τg
    return (DINobsmean .- DIN) ./ τg
end
function georesJac(p::Para)
    τg = p.τg
    return -I / τg
end

# Uptake
relu(x) = (x .≥ 0) .* x
drelu(x) = (x .≥ 0) .* 1.0
# Michaelis-Menten
mm(x, μ, k)     =  μ * x ./ (x .+ k)
∂mm_∂x(x, μ, k) =  μ * k ./ (x .+ k).^2
∂mm_∂μ(x, μ, k) =      x ./ (x .+ k)
∂mm_∂k(x, μ, k) = -μ * x ./ (x .+ k).^2
function uptake(DIN, p::Para)
    umax, ku = p.umax, p.ku
    return maskEup .* mm(relu(DIN), umax, ku)
end
# overload for numJac
function uptake(DIN::Array{<:Number,2}, p::Para)
    umax, ku = p.umax, p.ku
    return d₀(maskEup) * mm(relu(DIN), umax, ku)
end
# Uptake derivatives
function uptakeJac(DIN, p::Para)
    umax, ku = p.umax, p.ku
    return d₀(maskEup .* ∂mm_∂x(relu(DIN), umax, ku) .* drelu(DIN))
end
function ∂uptake_∂umax(DIN, p::Para)
    umax, ku = p.umax, p.ku
    return maskEup .* ∂mm_∂μ(relu(DIN), umax, ku)
end
function ∂uptake_∂ku(DIN, p::Para)
    umax, ku = p.umax, p.ku
    return maskEup .* ∂mm_∂k(relu(DIN), umax, ku)
end

# Remineralization
function remineralization(POM, p::Para)
    κ = p.κ
    return κ * POM
end
# Remineralization derivatives
function remineralizationJac(POM, p::Para)
    κ = p.κ
    return κ * I
end
function ∂remineralization_∂κ(POM, p::Para)
    κ = p.κ
    return POM
end

# Indices for DIN and POM (needed for Rate of change F(x,p))
const iDIN = 1:nwet
const iPOM = iDIN .+ nwet
function unpackx(x)
    DIN = x[iDIN]
    POM = x[iPOM]
    return DIN, POM
end
# add method to deal with arrays for numJac
function unpackx(x::Array{<:Number,2})
    DIN = x[iDIN,:]
    POM = x[iPOM,:]
    return DIN, POM
end

# PFD transport (needed for Rate of change F(x,p))
const S1 = buildPFD(ones(nwet), DIV, Iabove)
const Sz = buildPFD(ztop, DIV, Iabove)
function S(p::Para)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S1 + w′ * Sz
end

# Rate of change F(x,p)
function F(x::Vector{Tx}, p::Para{Tp}) where {Tx, Tp}
    DIN, POM = unpackx(x)
    u = uptake(DIN, p)
    r = remineralization(POM, p)
    foo = zeros(promote_type(Tx, Tp), size(x))
    foo[iDIN] .=  -T   * DIN - u + r + geores(DIN, p)
    foo[iPOM] .= -S(p) * POM + u - r
    return foo
end
# For numJac
function F(x::Array{<:Number,2}, p::Para)
    DIN, POM = unpackx(x)
    u = uptake(DIN, p)
    r = remineralization(POM, p)
    foo = zeros(eltype(x), size(x))
    foo[iDIN, :] .=  -T   * DIN - u + r + geores(DIN, p)
    foo[iPOM, :] .= -S(p) * POM + u - r
    return foo
end

# Jacobian of f with respect to x
function ∇ₓF(x, p::Para)
    DIN, POM = unpackx(x)
    uJac = uptakeJac(DIN, p)
    rJac = remineralizationJac(POM, p)
    foo = [ -T   - uJac + georesJac(p)    rJac       ;
                  +uJac                  -rJac - S(p)]
    dropzeros!(foo)
    return foo
end

using Flatten
"""
    ∇ₚF(x, p::Para)

Evaluates the jacobian of `f` with respect to `p`.
Concatenates `∇ₚF(x, p::Para, s::Symbol)` for all optimizable parameter symbols `s`.
"""
∇ₚF(x, p::Para) = hcat((∇ₚF(x, p, s) for s in fieldnameflatten(p))...)

"""
    ∇ₚF(x, p::Para, s::Symbol)

Evaluates the derivative of `f` with respect to `p.s` where `s` is the name of the parameter (of type `Symbol`).

You should fill this function with all the first derivatives of f with respoect to each parameter.
Called without the symbol `s`, this function will loop through all the optimizable parameters and create the corresponding Jacobian matrix.
"""
function ∇ₚF(x, p::Para, s::Symbol)
    DIN, POM = unpackx(x)
    foo = zeros(promote_type(eltype(x), eltype(p)), 2nwet)
    if s == :umax
        foo[iPOM] .= ∂uptake_∂umax(DIN, p)
        foo[iDIN] .= -foo[iPOM]
    elseif s == :ku
        foo[iPOM] .= ∂uptake_∂ku(DIN, p)
        foo[iDIN] .= -foo[iPOM]
    elseif s == :w₀
        foo[iPOM] .= -S1 * POM
    elseif s == :w′
        foo[iPOM] .= -Sz * POM
    elseif s == :κ
        foo[iDIN] .= ∂remineralization_∂κ(POM, p)
        foo[iPOM] .= -foo[iDIN]
    else
        error("There is no $s parameter")
    end
    return foo
end

