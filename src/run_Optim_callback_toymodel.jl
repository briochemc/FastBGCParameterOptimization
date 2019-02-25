# Main script to plot paper figures

# 1. Always load the packages
include("load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)

include("build_6BoxModel_circulation.jl")
using .SixBoxModel: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

# include("OCIM1.jl")
# using .OCIM1: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

# The above must define everything needed by the functions

# load biogeochmistry parameters
include("bgc_parameters.jl")

# load automatic differentiation stuff
include("AutomaticDifferentiation.jl")

# load biogeochmistry functions
include("bgc_functions.jl")

# load cost functions
include("cost_functions.jl")

# Run benchmark on q, Dq, D2q, etc.
include("save_Optim_callback_toymodel_data.jl")
