# Biogeochemistry parameters
module Parameters

# Packages for parameters
using FieldDefaults, Flatten, FieldMetadata, Unitful, Distributions
import FieldDefaults: get_default
import FieldMetadata: @units, units, @prior, prior, @description, description
import Flatten: flattenable

import TransportMatrixTools: AbstractPara

# Define some metadata
@metadata printunits nothing
@metadata description ""
@metadata latexSymbol ""
# Define the year unit (not in Unitful)
@unit(yr, "yr", Year, 365.0u"d", true)
Unitful.register(@__MODULE__)

const spd = 24*60*60
"""
    Para{U} <: AbstractPara{U}

Biogeochemical parameters (type).
You can create parameters `p` by using the default constructor `p = Para()`, which will have the default values.
Or you can specify some of the parameter values via `p = Para(p₁ = <value1>, p₂ = <value2>, ...)`.
Make sure you know what the parameters are!
Each parameter in `p` comes with a bunch of metadata for each field `f`:
- `latexSymbol(p, f)` is the LaTeX string of the unicode symbol (`String`)
- `description(p, f)` describes the parameter (`String`)
- `flattenable(p, f)` — is p.f optimizable? (`Boolean`)
- `prior(p, f)` gives the prior distribution (`Distribution type`)
- `printunits(p, f)` gives the unit used for `show(p)`
- `units(p, f)` gives the unit used in the model (SI)
- `default(p, f)` gives the default value (type `U`)
Modify this part of the code if you need new/different parameters!
"""
@latexSymbol @description @flattenable @prior @printunits @units @default_kw mutable struct Para{U} <: AbstractPara{U}
    umax::U |  1e-3 / spd | u"mmol/m^3/s" | u"mmol/m^3/d" | LN(1e-3 / spd, 1e-3 / spd) | true  | "Maximum uptake rate (Michaelis-Menten)"      | "\\mathbf{u}_\\mathrm{max}"
    ku  ::U |        1e-1 | u"mmol/m^3"   | u"mmol/m^3"   | LN(1e-1, 0.5e-1)           | true  | "Half-saturation constant (Michaelis-Menten)" | "k_\\mathbf{u}"
    w₀  ::U |     1 / spd | u"m/s"        | u"m/d"        | LN(1 / spd, 0.5 / spd)     | true  | "Sinking velocity at surface"                 | "w_0"
    w′  ::U |     1 / spd | u"s^-1"       | u"d^-1"       | LN(1 / spd, 0.5 / spd)     | true  | "Vertical gradient of sinking velocity"       | "w'"
    κ   ::U |  0.25 / spd | u"s^-1"       | u"d^-1"       | LN(0.25 / spd, 0.1 / spd)  | true  | "Remineralization rate"                       | "\\kappa"
    τg  ::U | 365e6 * spd | u"s"          | u"yr"         | nothing                    | false | "Geological Restoring"                        | "\\tau_\\mathrm{geo}"
    ω   ::U |        1e-4 | u"1"          | u"1"          | nothing                    | false | "Relative weight of params in cost"           | "\\omega"
end

"""
    LN(m, s)

Gives the `LogNormal` distribution that has a mean `m` and standard deviation `s`.
(I use the standard deviation rather than the variance to avoid unit-conversion confusion.)
"""
LN(m, s) = LogNormal(μ_LogNormal(m, s), σ_LogNormal(m, s))
μ_LogNormal(m, s) = log(m / sqrt(1 + (s / m)^2))
σ_LogNormal(m, s) = sqrt(log(1 + (s / m)^2))


"""
    optimizable_parameters

Tuple (constant) of the symbols of the optimizable parameters.
Uses the `Flatten.jl` package.
"""
const optimizable_parameters = fieldnameflatten(Para())

"""
    m

Number of optimizable parameters (constant).
"""
const m = length(fieldnameflatten(Para()))

"""
    m_all

Number of parameters including non-optimizable ones (constant).
"""
const m_all = length(fieldnames(typeof(Para())))

Base.eltype(::Para{U}) where U = U

"""
    length(p::Para)

Returns the number of optimizable parameters of `p`.
"""
Base.length(p::Para) = length(fieldnameflatten(p))

# Defining an iterator for the parameters to be able to `collect` it into a vector
Base.iterate(p::Para, i=1) = i > m_all ? nothing : (getfield(p, i), i + 1)

# Convert p to a vector and vice versa
Base.vec(p::Para) = collect((p...,))
Para(v::Vector) = Para(v...)
Base.copy(p::Para) = Para(vec(p)...)

# read non-real part (for update of init)
using DualNumbers, HyperDualNumbers
DualNumbers.realpart(p::Para{Dual{Float64}}) = Para(DualNumbers.realpart.(vec(p)))
Base.real(p::Para{Complex{Float64}}) = Para(real.(vec(p)))
HyperDualNumbers.realpart(p::Para{Hyper{Float64}}) = Para(HyperDualNumbers.realpart.(vec(p)))

# Overloads for parameters
Base.:+(p₁::Para, p₂::Para) = Para(vec(p₁) .+ vec(p₂))
Base.:*(p₁::Para, p₂::Para) = Para(vec(p₁) .* vec(p₂))
Base.:*(s::Number, p::Para) = Para(s .* vec(p))
Base.:*(p::Para, s::Number) = Para(s .* vec(p))
Base.isapprox(p₁::Para, p₂::Para) = isapprox(vec(p₁), vec(p₂))
Base.:(==)(p₁::Para, p₂::Para) = vec(p₁) == vec(p₂)
strerror = "No! Can't access the parameters at this index!"
Base.getindex(p::Para, i::Int) = i < 1 || i > m ? error(strerror) : getfield(p, i)
Base.setindex!(p::Para, v, i) = i < 1 || i > m ? error(strerror) : setfield!(p::Para, i, v)

# Convert p to λ and vice versa, needed by TransportMatrixTools!
optvec(p::Para) = flatten(Vector, p)
Base.:+(p::Para, v::Vector) = p + Flatten.reconstruct(Para((zeros(eltype(v), m_all))...), v)

# For shifting initial state by ∇s*∇p when p is updated
Base.:-(p₁::Para, p₂::Para) = Flatten.reconstruct(p₁, optvec(p₁) .- optvec(p₂))
Base.:*(∇s::Array, p::Para) = ∇s * optvec(p)

"""
    p₀

The (constant) default values of non-optimized parameters.
"""
const p₀ = Para(ω = 1e-4)
str_out = "_default"

"""
    μobs

The (constant) mean of the log of the observed parameters (the μ of the lognormal prior).
"""
const μobs = [meanlogx.(metaflatten(p₀, prior))...]

"""
    σ²obs

The (constant) variance of the log of the observed parameters (the σ² of the lognormal prior).
"""
const σ²obs = [varlogx.(metaflatten(p₀, prior))...]

"""
    p2λ(p::Para)

Returns the `λ` that corresponds to `p`.
`p2λ` and `λ2p` allow for functionality:
E.g., conditions like being positive can be imposed using `exp`, etc.
"""
p2λ(p::Para) = log.(optvec(p)) - μobs
const λ₀ = p2λ(p₀)
∇p2λ(p::Para) = optvec(p).^(-1)
λ2p(λ) = optPara(exp.(λ + μobs))
∇λ2p(λ) = exp.(λ + μobs)
∇²λ2p(λ) = exp.(λ + μobs)

function optPara(v::Vector{U}) where U
    if U == eltype(p₀)
        return Flatten.reconstruct(p₀, v)
    else
        p = Para(convert(Vector{U}, vec(p₀)))
        return Flatten.reconstruct(p, v)
    end
end

#=============================================
    Functions for printing
=============================================#

using Printf
import Base: show
function Base.show(io::IO, p::Para)
    println(typeof(p))
    for f in fieldnames(typeof(p))
        v = getfield(p, f)
        punit = printunits(p, f)
        if punit == 1
            val, ppunit = v, ""
        else
            val, ppunit = ustrip(uconvert(punit, v * units(p, f))), unicodify(punit)
        end
        @printf io "%6s = %8.2e %s\n" f val ppunit
    end
end
Base.show(io::IO, ::MIME"text/plain", p::Para) = Base.show(io, p)

function unicodify(U::Unitful.Units)
    str = string(U)
    str = replace(str, r"\^-1" => s"⁻¹")
    str = replace(str, r"\^-2" => s"⁻²")
    str = replace(str, r"\^-3" => s"⁻³")
    str = replace(str, r"\^1" => s"¹")
    str = replace(str, r"\^2" => s"²")
    str = replace(str, r"\^3" => s"³")
    return str
end

#=============================================
    Functions for printing to LaTeX table TODO
=============================================#

function print_LaTeX_table(p::Para)
    println("Latex table for the parameters:")
    println("\\begin{table*}[t]")
    println("    \\centering")
    println("    \\begin{tabular}{lllll}")
    println("        Parameter & Description & Value & Unit & Optimized? \\\\ \\hline")
    for f in fieldnames(typeof(p))
        LaTeX_line = "        " # Indent
        LaTeX_line = LaTeX_line * "\$" * latexSymbol(p, f) * "\$" # Parameter Symbol
        LaTeX_line = LaTeX_line * " & " * description(p, f) # Parameter Description
        val = uconvert(printunits(p, f), getfield(p, f) * units(p, f)) # Value in print units
        LaTeX_line = LaTeX_line * " & \$" * latexify(ustrip(val)) * "\$"
        LaTeX_line = LaTeX_line * " & " * latexify(printunits(p, f)) # print units
        LaTeX_line = LaTeX_line * " & " * latexify(flattenable(p, f)) * "\\\\"
        println(LaTeX_line)
    end
    println("        \\label{T:Parameters}")
    println("    \\end{tabular}")
    println("    \\caption{")
    println("    }")
    println("\\end{table*}")
end

function latexify(U::Unitful.Units)
    str = string(U)
    str = replace(str, r"\^-(?<exp>\d+)" => s"$^{-\g<exp>}$") # add brackets around exponents
    str = replace(str, r"dy" => s"d") # Replace "dy" with "d" for days
    str = replace(str, r"\s" => s"\\,") # replace spaces with small spaces
    return str
end

latexify(optimized::Bool) = optimized ? "yes" : "no"

function latexify(val)
    str = @sprintf("%.2g", val)
    str = replace(str, r"e\+0+(?<exp>\d+)" => s" \\times 10^{\g<exp>}")
    str = replace(str, r"e\-0+(?<exp>\d+)" => s" \\times 10^{-\g<exp>}")
    return str
end

# TODO use DataFrames, but solve export conflicts with, e.g., `compress`
#using DataFrames
#function tablify(p::Para)
#    s = collect(latexSymbol(p)) # Parameter Symbol
#    d = collect(description(p)) # Parameter Description
#    v = [ustrip(uconvert(printunits(p, f), getfield(p, f) * units(p, f))) for f in fieldnames(typeof(p))] # Value in print units
#    u = collect(printunits(p)) # print units
#    o = collect(flattenable(p))
#    return DataFrame(Parameter=s, Description=d, Value=v, Unit=u, optimized=o)
#end

# # macro for printing p::Para?
# function printPara(p::Para)
#     for param in p
#         Expr(:call, :println, string(param), " = ", esc(param))
#     end
# end


end
