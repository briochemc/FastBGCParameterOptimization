#!/bin/bash
#$ -N run_methods_benchmark
#$ -q pub64
#$ -pe openmp 8-64
#$ -l mem=64G
#$ -m beas

# Load the julia module
module load julia/1.1.0

# Go to the root folder on greenplanet
cd /data/users/pasquieb/Projects/FastBGCParameterOptimization

# Set DataDeps environment variable to download without asking
DATADEPS_ALWAYS_ACCEPT = true

# Run it!
julia src/run_methods_benchmark.jl
