# Load circulation
const wet3d, grd, T_Circulation = Circulation.load()

# Define useful constants and arrays
const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)
# And matrices
const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet)

# Generate biogeochmistry Parameters
include("bgc_parameters.jl")
# Define biogeochmistry functions
include("bgc_functions.jl")
# Define mismatch functions
include("cost_functions.jl")
# Load methods
include("load_methods.jl")
# Use methods to define objective, gradient, and Hessian
include("use_methods.jl")
# Set extra constants like τstop, x₀, etc. (Maybe should be elsewhere)
include("meta_constants.jl")

