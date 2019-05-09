# script to load all the packages needed

# activate current Julia project
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# use TransportMatrixTools package that I dev
using AIBECS

# Matrix stuff
using SparseArrays, SuiteSparse, LinearAlgebra
SuiteSparse.UMFPACK.umf_ctrl[8] = 0 # turn iterative refinements off

# # pretty print
# using Printf
# 
# # Dual and Hyperdual numbers (and my tools)
# using DualNumbers, HyperDualNumbers
# using DualMatrixTools, HyperDualMatrixTools

# Loading and saving data
using JLD2
# using JLD2, MAT
# Note: I put MAT only where it is called because problems compiling MAT.jl in Julia v1.0+

# Optimization package
using Optim

# Package to load and grid WOA data
# using WorldOceanAtlasTools (should be part of AIBECS now)

# Benchmarking
using BenchmarkTools

