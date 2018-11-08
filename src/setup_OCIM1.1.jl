const spd = 24 * 60 * 60.0 # Think about using Unitful.jl

function load_OCIM_mat_file()
    println("loading OCIM with MAT")
    vars = matread(homedir() * "/Dropbox/OCIM1.1/CTL.mat")
    T = -365spd * vars["output"]["TR"]
    grd = vars["output"]["grid"]
    wet3d = vars["output"]["M3d"]
    return wet3d, grd, T
end


if isfile(joinpath(pwd(), "data", "OCIM1.1_CTL.jld2"))
    println("loading OCIM with JLD2")
    @load joinpath(pwd(), "data", "OCIM1.1_CTL.jld2") wet3d grd T
else
    wet3d, grd, T = load_OCIM_mat_file()
    println("saving OCIM with JLD2")
    @save joinpath(pwd(), "data", "OCIM1.1_CTL.jld2") wet3d grd T
end


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

# load observations - Recode the fetching and regridding in Julia!
@load (homedir() * "/.julia/v0.6/OCIMtools/data/OCIM1.1_SiObs.jld2") xobs
const DSiobs = xobs
const DSimean = vmean(DSiobs)

# Required for bgc functions
const maskEup = z .< 85.0 # Euphotic zone definition


