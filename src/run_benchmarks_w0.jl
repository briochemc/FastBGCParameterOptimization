# Main script to plot paper figures

# 1. Always load the packages
include("load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)

# include("build_6BoxModel_circulation.jl")
# using .SixBoxModel: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

include("OCIM1.jl")
using .OCIM1: T, wet3d, grd, spd, nwet, DIV, Iabove, ztop, DINobs , vnorm², maskEup, DINobsmean, Dvnorm²

# The above must define everything needed by the functions

# load biogeochmistry parameters
include("bgc_parameters.jl")

# load automatic differentiation stuff
include("AutomaticDifferentiation.jl")

# load biogeochmistry functions
include("bgc_functions.jl")

# load cost functions
include("cost_functions.jl")

# set meta constants: p₀, λ₀, τstop, x₀, etc.
include("meta_constants_w0.jl")

# Run benchmark on q, Dq, D2q, etc.
include("save_D2q_benchmark_data.jl")
