# Biogeochemistry parameters

# BGC parameters
# Define some metadata
@metadata printunits nothing
@metadata description ""
@metadata latexSymbol ""
# Define the year unit (not in Unitful)
@unit(yr, "yr", Year, 365.0u"d", true)
Unitful.register(@__MODULE__)
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
@latexSymbol @description @flattenable @prior @printunits @units @default_kw struct Para{U} <: AbstractPara{U}
    τu::U | 500.0 * spd | u"s"    | u"d"    | LN(250 * spd, 100 * spd)  | true  | "Specific uptake rate timescale"        | "\\tau_\\mathbf{u}"
    w₀::U |     1 / spd | u"m/s"  | u"m/d"  | LN(1 / spd, 0.5 / spd)    | true  | "Sinking velocity at surface"           | "w_0"
    w′::U |     1 / spd | u"s^-1" | u"d^-1" | LN(1 / spd, 0.5 / spd)    | true  | "Vertical gradient of sinking velocity" | "w'"
    κ::U  |  0.25 / spd | u"s^-1" | u"d^-1" | LN(0.25 / spd, 0.1 / spd) | true  | "Remineralization rate"                 | "\\kappa"
    τg::U | 365e6 * spd | u"s"    | u"yr"   | nothing                   | false | "Geological Restoring"                  | "\\tau_\\mathrm{geo}"
end

"""
    LN(m, s)

Gives the `LogNormal` distribution that has a mean `m` and standard deviation `s`.
(I use the standard deviation rather than the variance to avoid unit-conversion confusion.)
"""
LN(m, s) = LogNormal(μ_LogNormal(m, s), σ_LogNormal(m, s))

"""
    μ_LogNormal(m, s)

Gives the mean, `μ`, of the log of the lognormal distribution with mean `m` and standard deviation `s`.
"""
μ_LogNormal(m, s) = log(m / sqrt(1 + (s / m)^2))

"""
    σ_LogNormal(m, s)

Gives the standard deviation, `σ`, of the log of the lognormal distribution with mean `m` standard deviation `s`.
"""
σ_LogNormal(m, s) = sqrt(log(1 + (s / m)^2))

"""
    p₀

The (constant) default values of non-optimized parameters.
"""
const p₀ = Para()

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
    optimizable_parameters

Tuple (constant) of the symbols of the optimizable parameters.
Uses the `Flatten.jl` package.
"""
const optimizable_parameters = fieldnameflatten(p₀)

"""
    npopt

Number of optimizable parameters (constant).
"""
const npopt = length(optimizable_parameters)

"""
    np

Number of parameters including non-optimizable ones (constant).
"""
const np = length(fieldnames(typeof(p₀)))

"""
    eltype(::Para{U})

Returns the element type of the parameters.
"""
Base.eltype(::Para{U}) where U = U

"""
    length(p::Para)

Returns the number of optimizable parameters of `p`.
"""
Base.length(p::Para) = length(fieldnameflatten(p₀))

# Defining an iterator for the parameters to be able to `collect` it into a vector
Base.iterate(p::Para, i=1) = i > np ? nothing : (getfield(p, i), i + 1)

# Convert p to a vector and vice versa
Base.vec(p::Para) = collect((p...,))
Para(v::Vector) = Para(v...)

# read non-real part (for update of init)
DualNumbers.realpart(p::Para{Dual{Float64}}) = Para(DualNumbers.realpart.(vec(p)))
Base.real(p::Para{Complex{Float64}}) = Para(real.(vec(p)))
HyperDualNumbers.realpart(p::Para{Hyper{Float64}}) = Para(HyperDualNumbers.realpart.(vec(p)))

# Overload +, -, and * for parameters
Base.:+(p₁::Para, p₂::Para) = Para(vec(p₁) .+ vec(p₂))
Base.:-(p₁::Para, p₂::Para) = Para(vec(p₁) .- vec(p₂))
Base.:*(p₁::Para, p₂::Para) = Para(vec(p₁) .* vec(p₂))
Base.:*(s::Number, p::Para) = Para(s .* vec(p))
Base.:*(p::Para, s::Number) = Para(s .* vec(p))
Base.isapprox(p₁::Para, p₂::Para) = isapprox(vec(p₁), vec(p₂))

# Convert p to λ and vice versa, needed by TransportMatrixTools!
optvec(p::Para) = flatten(Vector, p)

"""
    p2λ(p::Para)

Returns the `λ` that corresponds to `p`.
`p2λ` and `λ2p` allow for functionality:
E.g., conditions like being positive can be imposed using `exp`, etc.
"""
p2λ(p::Para) = log.(optvec(p)) - μobs


"""
    λ₀

Default `λ` corresponding to `p₀` (constant).
"""
const λ₀ = p2λ(p₀)

"""
    Dp2λ(p::Para)

Returns the gradient of `λ` with respect to `p`.
"""
Dp2λ(p::Para) = optvec(p).^(-1)

function optPara(v::Vector{U}) where U
    if U == eltype(p₀)
        return Flatten.reconstruct(p₀, v)
    else
        p = Para(convert(Vector{U}, vec(p₀)))
        return Flatten.reconstruct(p, v)
    end
end

"""
    λ2p(λ)

Returns the `p` that corresponds to `λ`.
See `p2λ`.
"""
λ2p(λ) = optPara(exp.(λ .+ μobs))

"""
    Dλ2p(λ)

Returns the gradient of `p` with respect to `λ`.
"""
Dλ2p(λ) = exp.(λ .+ μobs)

"""
    D2λ2p(λ)

Returns the Hessian of `p` with respect to `λ`.
"""
D2λ2p(λ) = exp.(λ .+ μobs)

"""
    show(io::IO, p::Para)

shows the parameters after applying the base unit and converting to printing unit.

`show(p)` will show all the parameters.

`show(IOContext(stdout, :compact => true), p)` will only show the optimizable parameters.
"""
function Base.show(io::IO, p::Para; preprint="")
    println(preprint * "Parameter values:")
    compact = get(io, :compact, false)
    for f in fieldnames(typeof(p))
        val = uconvert(printunits(p, f), getfield(p, f) * units(p, f))
        (~compact || flattenable(p, f)) && println(preprint * "│ $f = $val")
    end
end

"""
    print_LaTeX_table(p::Para)

prints a LaTeX table of the parameters applying the printing unit.
"""
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

"""
    latexify(val)

Returns the LaTeX string for the value `val` to be used in the LaTeX table of parameters.
See also: [`print_LaTeX_table`](@ref)
"""
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

