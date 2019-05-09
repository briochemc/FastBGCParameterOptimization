# Load the packages
include("load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)
# Circulation = OCIM1
Circulation = SixBoxModel

const wet3d, grd, T_Circulation = Circulation.load()


const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)

const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet)
#=
spd, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²?
=#


# load biogeochmistry parameters
include("bgc_parameters.jl")
# load biogeochmistry functions
include("bgc_functions.jl")
# load cost functions
include("cost_functions.jl")
#=
Replace by AIBECS style
except ∇ₓF, ∇ₓf, ∇ₚF, and ∇ₚf must be analytical for DUAL, CSD, and HYPER methods
=#

# load methods
include("load_methods.jl")

# load derivatives functions
include("use_methods.jl")

# set meta constants: p₀, λ₀, τstop, x₀, etc.
include("meta_constants.jl")

