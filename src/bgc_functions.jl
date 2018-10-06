# Biogeochemistry functions

# BGC parameters
# Define some metadata
@metadata printunits nothing
@metadata describe ""
# Define the year unit (not in Unitful)
@unit(yr, "yr", Year, 365.0u"d", true)
Unitful.register(@__MODULE__)
"""
    Para

Biogeochemical parameters (type).
You can create parameters `p` by using the default constructor `p = Para()`, which will have the default values.
Or you can specify some of the parameter values via `p = Para(p₁ = <value1>, p₂ = <value2>, ...)`.
Make sure you know what the parameters are!
Each parameter in `p` comes with a bunch of metadata for each field `f`:
- `flattenable(p, f)` — is p.f optimizable? (boolean)
- `describe(p, f)` describes the parameter (string)
- `printunits(p, f)` gives the unit used for `show(p)`
- `units(p, f)` gives the unit used in the model (SI)
- `default(p, f)` gives the default value
Modify this part of the code if you need new/different parameters!
"""
@describe @flattenable @printunits @units @default_kw struct Para{U}
    τu::U |  50.0 * spd | u"s"    | u"d"    | true  | "Specific uptake rate timescale"
    w₀::U |     1 / spd | u"m/s"  | u"m/d"  | true  | "Sinking velocity at surface"
    w′::U |     1 / spd | u"s^-1" | u"d^-1" | false | "Vertical gradient of sinking velocity"
    κ::U  |  0.25 / spd | u"s^-1" | u"d^-1" | false | "Remineralization rate"
    τg::U | 365e6 * spd | u"s"    | u"yr"   | false | "Geological Restoring"
end
const p₀ = Para() # p₀ will hold the default values of non-optimized parameters
const pobs = Para(τu = 50.0 * spd, w₀ = 100 / spd)
const optimizable_parameters = fieldnameflatten(p₀)
const npopt = length(optimizable_parameters)
const all_parameters = fieldnames(typeof(p₀))
const np = length(all_parameters)
# Overload eltype to figure out the type of the parameters
# Required because now parameters are in a struct :)
Base.eltype(::Para{U}) where U = U
Base.length(p::Para) = length(fieldnameflatten(p₀)) # This is dangerous! But it may work :)
# Defining an iterator for the parameters
Base.iterate(p::Para, i=1) = i > np ? nothing : (getfield(p, i), i + 1)
# Convert p to a vector and vice versa
Base.vec(p::Para) = collect((p...,))
Para(v::Vector) = Para(v...)
# Overload +, -, and * for parameters
Base.:+(p₁::Para, p₂::Para) = Para(vec(p₁) .+ vec(p₂))
Base.:-(p₁::Para, p₂::Para) = Para(vec(p₁) .- vec(p₂))
Base.:*(p₁::Para, p₂::Para) = Para(vec(p₁) .* vec(p₂))
Base.:*(s::Number, p::Para) = Para(s .* vec(p))
# Convert p to λ and vice versa, needed by TransportMatrixTools!
optvec(p::Para) = flatten(Vector, p)
p2λ(p::Para) = log.(optvec(p) ./ optvec(pobs))
Dp2λ(p::Para) = optvec(p).^(-1)
function optPara(v::Vector{U}) where U
    if U == eltype(p₀)
        return Flatten.reconstruct(p₀, v)
    else
        p = Para(convert(Vector{U}, vec(p₀)))
        return Flatten.reconstruct(p, v)
    end
end
λ2p(λ) = optPara(exp.(λ) .* optvec(pobs))
Dλ2p(λ) = exp.(λ) .* optvec(pobs)
D2λ2p(λ) = exp.(λ) .* optvec(pobs)

"""
    show(io::IO, p::Para)

shows the parameters after applying the base unit and converting to printing unit.

`show(p)` will show all the parameters.

`show(IOContext(stdout, :compact => true), p)` will only show the optimizable parameters.
"""
function Base.show(io::IO, p::Para)
    println("Parameter values:")
    compact = get(io, :compact, false)
    for f in fieldnames(typeof(p))
        val = uconvert(printunits(p, f), getfield(p, f) * units(p, f))
        (~compact || flattenable(p, f)) && println("| $f = $val")
    end
end


# # macro for printing p::Para?
# function printPara(p::Para)
#     for param in p
#         Expr(:call, :println, string(param), " = ", esc(param))
#     end
# end


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
# Remineralization derivative
function remineralizationJac(PSi, p::Para)
    κ = p.κ
    return κ * I
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

# # Derivative of f with respect to p
# function Dpf(x, p)
#     DSi, PSi = unpackx(x)
#     foo = zeros(promote_type(eltype(x), eltype(p)), 2nwet, np)
#     foo[iPSi, 1] .= DuptakeDτu(DSi, p)
#     foo[iDSi, 1] .= -foo[iPSi, 1]
#     foo[iPSi, 2] .= -S1 * PSi
#     foo[iPSi, 3] .= -Sz * PSi
#     return foo

# end



Dpf(x, p::Para) = hcat((Dpf(x, p, popt) for popt in optimizable_parameters)...)

"""
    Dpf(x, p::Para, s::Symbol)

Evaluates the derivative of `s` with respect to `p.s` where `s` is the name of the parameter (of type `Symbol`).

You should fill this function with all the first derivatives of f with respoect to each parameter.
Called without the symbol `s`, this function will loop through all the optimizable parameters and create the corresponding matrix.
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
        foo[iPSi] .= -Sz * PSi
    else
        error("There is no $s parameter")
    end
    return foo
end

