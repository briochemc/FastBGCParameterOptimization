#!/bin/bash

#SBATCH --job-name=run_methods_benchmark
#SBATCH --output=cluster_output/run_methods_benchmark.out
#SBATCH --error=cluster_output/run_methods_benchmark.err
#SBATCH --partition=brd2.4,has2.5,ilg2.3,m-c1.9,m-c2.2,nes2.8,sib2.9
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64GB

# Load the julia module
ml purge
ml julia/1.1.0

# Go to the root folder on greenplanet
cd /DFS-L/DATA/moore/pasquieb/Projects/FastBGCParameterOptimization

# Set DataDeps environment variable to download without asking
DATADEPS_ALWAYS_ACCEPT = true

# Run it!
julia src/run_methods_benchmark.jl
