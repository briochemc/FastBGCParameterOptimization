# These are mostly commands to copy paste for now as I am figuring this all out...

# Install git
# sudo apt-get install git

# Install Julia
# wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz
# tar xvzf julia-1.1.0-linux-x86_64.tar.gz

# Git clone my project
# git clone https://github.com/briochemc/FastBGCParameterOptimization.git

# cd into project
# cd FastBGCParameterOptimization

# Run it
# /home/briochemc/julia-1.1.0/bin/julia src/run_methods_benchmark.jl

# Build NetCDF if it bugs
# /home/briochemc/julia-1.1.0/bin/julia -e 'using Pkg; Pkg.add("NetCDF"); Pkg.build("NetCDF")'

# Go to the root folder on Google Cloud

# Set DataDeps environment variable to download without asking
DATADEPS_ALWAYS_ACCEPT = true

# Run it!
julia src/run_methods_benchmark.jl

