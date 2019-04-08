# Load the packages
include("load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)

include("build_6BoxModel_circulation.jl")
Circulation = SixBoxModel
# include("OCIM1.jl")
# Circulation = OCIM1
using .Circulation: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

# load biogeochmistry parameters
include("bgc_parameters.jl")
using .Parameters: Para, p₀, λ₀, p2λ, ∇p2λ, ∇²λ2p, ∇λ2p, λ2p

# load biogeochmistry functions
include("bgc_functions.jl")

# load cost functions
include("cost_functions.jl")

# load methods
include("load_methods.jl")

# load derivatives functions
include("use_methods.jl")

# set meta constants: p₀, λ₀, τstop, x₀, etc.
include("meta_constants.jl")

