# Biogeochemistry parameters

# BGC parameters
# Define some metadata
@metadata printunits nothing
@metadata describe ""
@metadata latexSymbol ""
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
@latexSymbol @describe @flattenable @printunits @units @default_kw struct Para{U}
    τu::U | 500.0 * spd | u"s"    | u"d"    | true  | "Specific uptake rate timescale"        | "\\tau_\\mathbf{u}"
    w₀::U |     1 / spd | u"m/s"  | u"m/d"  | true  | "Sinking velocity at surface"           | "w_0"
    w′::U |     1 / spd | u"s^-1" | u"d^-1" | false | "Vertical gradient of sinking velocity" | "w'"
    κ::U  |  0.25 / spd | u"s^-1" | u"d^-1" | false | "Remineralization rate"                 | "\\kappa"
    τg::U | 365e6 * spd | u"s"    | u"yr"   | false | "Geological Restoring"                  | "\\tau_\\mathrm{geo}"
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
# read non-real part (for update of xinit)
realpart(p::Para) = Para(realpart.(vec(p)))
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
        (~compact || flattenable(p, f)) && println("│ $f = $val")
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
        LaTeX_line = LaTeX_line * " & " * describe(p, f) # Parameter Description
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

function latexify(val::R) where R<:Real
    str = @sprintf("%.2g", val)
    str = replace(str, r"e\+0+(?<exp>\d+)" => s" \\times 10^{\g<exp>}")
    str = replace(str, r"e\-0+(?<exp>\d+)" => s" \\times 10^{-\g<exp>}")
    return str
end

using DataFrames
function tablify(p::Para)
    s = collect(latexSymbol(p)) # Parameter Symbol
    d = collect(describe(p)) # Parameter Description
    v = [ustrip(uconvert(printunits(p, f), getfield(p, f) * units(p, f))) for f in fieldnames(typeof(p))] # Value in print units
    u = collect(printunits(p)) # print units
    o = collect(flattenable(p))
    return DataFrame(Parameter=s, Description=d, Value=v, Unit=u, optimized=o)
end

# # macro for printing p::Para?
# function printPara(p::Para)
#     for param in p
#         Expr(:call, :println, string(param), " = ", esc(param))
#     end
# end

