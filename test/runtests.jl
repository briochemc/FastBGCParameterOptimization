using Test

# 1. Always load the packages
include("../src/load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)
# include("setup_OCIM1.1.jl")
include("../src/setup_4BoxModel.jl")
# The above must define everything needed by the functions

# load biogeochmistry functions
include("../src/bgc_functions.jl")

# load cost functions
include("../src/cost_functions.jl")

# run tests
include("test_Derivatives.jl")
