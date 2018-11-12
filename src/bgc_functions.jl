
# Geological restoring
function geores(DSi, p::Para)
    τg = p.τg
    return (DSimean .- DSi) ./ τg
end
function georesJac(p::Para)
    τg = p.τg
    return -I / τg
end

# Uptake
relu(x) = (x .≥ 0) .* x
Drelu(x) = (x .≥ 0) .* 1.0
function uptake(DSi, p::Para)
    τu = p.τu
    return maskEup .* relu(DSi - DSiobs) ./ τu
end
# overload for numJac
function uptake(DSi::Array{<:Number,2}, p::Para)
    τu = p.τu
    DSiobs2 = repeat(DSiobs, 1, size(DSi, 2))
    return d₀(maskEup) * relu(DSi - DSiobs2) ./ τu
end
# Uptake derivatives
function uptakeJac(DSi, p::Para)
    τu = p.τu
    return d₀(maskEup .* Drelu(DSi - DSiobs) ./ τu)
end
function Duptake_Dτu(DSi, p::Para)
    τu = p.τu
    return -uptake(DSi, p) ./ τu
end

# Remineralization
function remineralization(PSi, p::Para)
    κ = p.κ
    return κ * PSi
end
# Remineralization derivatives
function remineralizationJac(PSi, p::Para)
    κ = p.κ
    return κ * I
end
function Dremineralization_Dκ(PSi, p::Para)
    κ = p.κ
    return PSi
end

# Indices for DSi and PSi (needed for Rate of change f(x,p))
const iDSi = 1:nwet
const iPSi = iDSi .+ nwet
function unpackx(x)
    DSi = x[iDSi]
    PSi = x[iPSi]
    return DSi, PSi
end
# add method to deal with arrays for numJac
function unpackx(x::Array{<:Number,2})
    DSi = x[iDSi,:]
    PSi = x[iPSi,:]
    return DSi, PSi
end

# PFD transport (needed for Rate of change f(x,p))
const S1 = buildPFD(ones(nwet), DIV, Iabove)
const Sz = buildPFD(ztop, DIV, Iabove)
function S(p::Para)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S1 + w′ * Sz
end

# Rate of change f(x,p)
function f(x::Vector{Tx}, p::Para{Tp}) where {Tx, Tp}
    DSi, PSi = unpackx(x)
    u = uptake(DSi, p)
    r = remineralization(PSi, p)
    foo = zeros(promote_type(Tx, Tp), size(x))
    foo[iDSi] .=  -T   * DSi - u + r + geores(DSi, p)
    foo[iPSi] .= -S(p) * PSi + u - r
    return foo
end
# For numJac
function f(x::Array{<:Number,2}, p::Para)
    DSi, PSi = unpackx(x)
    u = uptake(DSi, p)
    r = remineralization(PSi, p)
    foo = zeros(eltype(x), size(x))
    foo[iDSi, :] .=  -T   * DSi - u + r + geores(DSi, p)
    foo[iPSi, :] .= -S(p) * PSi + u - r
    return foo
end

# Jacobian of f with respect to x
function fJac(x, p::Para)
    DSi, PSi = unpackx(x)
    uJac = uptakeJac(DSi, p)
    rJac = remineralizationJac(PSi, p)
    foo = [ -T   - uJac + georesJac(p)    rJac       ;
                  +uJac                  -rJac - S(p)]
    dropzeros!(foo)
    return foo
end

"""
    Dpf(x, p::Para)

Evaluates the jacobian of `f` with respect to `p`.
Concatenates `Dpf(x, p::Para, s::Symbol)` for all optimizable parameter symbols `s`.
"""
Dpf(x, p::Para) = hcat((Dpf(x, p, popt) for popt in optimizable_parameters)...)

"""
    Dpf(x, p::Para, s::Symbol)

Evaluates the derivative of `f` with respect to `p.s` where `s` is the name of the parameter (of type `Symbol`).

You should fill this function with all the first derivatives of f with respoect to each parameter.
Called without the symbol `s`, this function will loop through all the optimizable parameters and create the corresponding Jacobian matrix.
"""
function Dpf(x, p::Para, s::Symbol)
    DSi, PSi = unpackx(x)
    foo = zeros(promote_type(eltype(x), eltype(p)), 2nwet)
    if s == :τu
        foo[iPSi] .= Duptake_Dτu(DSi, p)
        foo[iDSi] .= -foo[iPSi]
    elseif s == :w₀
        foo[iPSi] .= -S1 * PSi
    elseif s == :w′
        foo[iPSi] .= -Sz * PSi
    elseif s == :κ
        foo[iDSi] .= Dremineralization_Dκ(PSi, p)
        foo[iPSi] .= -foo[iDSi]
    else
        error("There is no $s parameter")
    end
    return foo
end

