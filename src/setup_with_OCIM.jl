# Load the packages
include("load_packages.jl")

const wet3d, grd, T_OCIM = TransportMatrixTools.OCIM1.load()

# load biogeochmistry parameters
include("bgc_parameters.jl")
using .Parameters: Para, p₀, λ₀, p2λ, ∇p2λ, ∇²λ2p, ∇λ2p, λ2p, σ²obs, getindex, m, str_out

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

