#!/bin/bash
#$ -N run_Cassette_profile
#$ -q pub64
#$ -pe openmp 8-64
#$ -l mem=64G
#$ -m beas

# Load the julia module
module load julia/1.1.0

# Go to the root folder on UCI HPC
cd /data/users/pasquieb/Projects/FastBGCParameterOptimization

# Set DataDeps environment variable to download without asking
DATADEPS_ALWAYS_ACCEPT = true

# Run it!
julia src/run_with_Cassette_profile.jl
