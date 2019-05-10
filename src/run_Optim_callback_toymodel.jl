# Load the packages
include("load_packages.jl")

# Define the transport matrix for the circulation
Circulation = SixBoxModel # from AIBECS

include("setup.jl")

include("save_Optim_callback_data.jl")
