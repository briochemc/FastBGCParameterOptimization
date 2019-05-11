
include("load_packages.jl")

# Define the transport matrix for the circulation
Circulation = OCIM # from AIBECS

include("setup.jl")

include("save_Optim_callback_data.jl")
