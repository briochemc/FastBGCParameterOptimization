using Test

include("../src/load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)

include("../src/build_6BoxModel_circulation.jl")
Circulation = SixBoxModel

#include("../src/OCIM1.jl")
#Circulation = OCIM1

using .Circulation: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

# load biogeochmistry parameters
include("../src/bgc_parameters.jl")

# load automatic differentiation stuff
include("../src/AutomaticDifferentiation.jl")

# load biogeochmistry functions
include("../src/bgc_functions.jl")

# load cost functions
include("../src/cost_functions.jl")

# set meta constants: p₀, λ₀, τstop, x₀, etc.
include("../src/meta_constants.jl")

# run tests
include("test_Derivatives.jl")


