module OCIM1
#=
This module serves to load the OCIM matrix and grid.
Those are loaded from the public, persistant URL in FigShare.
The Julia JLD2 format version of the OCIM was created with the code
in the GitHub repository https://github.com/briochemc/MatricesForJAMES.
The citation material is in here but I could make it availabe to future JAMES users.
Need to decide if JAMES is the name I want to keep too.
=#

using TransportMatrixTools: buildv3d, d₀, buildDIV, buildIabove
using SparseArrays, SuiteSparse
using WorldOceanAtlasTools
using DataDeps, JLD2, FileIO

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end

# Create registry entry for OCIM in JLD2 format
register(
    DataDep(
        "MatricesForJAMES",
        """
        References:
        - DeVries, T. (2014), The oceanic anthropogenic CO2 sink: Storage, air‐sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
        - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, https://doi.org/10.1175/JPO-D-10-05011.1
        """,
        "http://files.figshare.com/14330492/OCIM1_CTL.jld2",
        sha2_256,
        fetch_method = fallback_download
    )
)

# Load it
const spd = 24 * 60 * 60.0 # Think about using Unitful.jl
function load_OCIM1_jld2_file()
    println("loading OCIM1 with JLD2")
    jld2_file = @datadep_str string("MatricesForJAMES/", "OCIM1_CTL.jld2")
    @load jld2_file vars
    T = -vars["output"]["TR"] / 365spd
    grd = vars["output"]["grid"]
    wet3d = vars["output"]["M3d"]
    return wet3d, grd, T
end

const wet3d, grd, T = load_OCIM1_jld2_file()

# For standard use with OCIM
const iwet = (LinearIndices(wet3d))[findall(!iszero, wet3d)] # replaces find(wet3d) :(
const nwet = length(iwet)
const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet)
const v3d = buildv3d(grd)
const v = v3d[iwet]
const vtot = sum(v)
const V = d₀(v)
vnorm²(x) = x' * V * x # Volume-weighted norm squared
vnorm²(x::Vector{Complex{Float64}}) = transpose(x) * V * x # Volume-weighted norm squared
Dvnorm²(x) = (2v .* x)' # Volume-weighted norm squared
Dvnorm²(x::Vector{Complex{Float64}}) = transpose(2v .* x) # Volume-weighted norm squared
vnorm(x) = sqrt(vnorm²(x))      # volume-weighted norm
vmean(x) = v'x / vtot   # volume-weighted mean
vmean(x::Vector{Complex{Float64}}) = transpose(v) * x / vtot   # volume-weighted mean
const z = grd["ZT3d"][iwet]
const ztop = grd["ZW3d"][iwet]

# Vector of observations
DINobs3d = WOA13_bin_to_grid(grd, "PO4", "annual", "1°", "mean")
const DINobs = DINobs3d[iwet]
const DINobsmean = vmean(DINobs)

# Required for bgc functions
const maskEup = z .< 85.0 # Euphotic zone definition

export T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

end

