# script to load all the packages needed

# activate current Julia project
using Pkg
Pkg.activate(".")

# use TransportMatrixTools package that I dev
using TransportMatrixTools

# Matrix stuff
using SparseArrays, SuiteSparse, LinearAlgebra

# pretty print
using Printf

# Differentiating sutff
using Calculus
using DualNumbers, HyperDualNumbers
using DualMatrixTools, HyperDualMatrixTools

# Loading and saving data
using JLD2
# using JLD2, MAT
# Note: I put MAT only where it is called because problems compiling MAT.jl in Julia v1.0+

# Optimization package
using Optim

# Packages for parameters
using FieldDefaults, Flatten, FieldMetadata, Unitful, Distributions
import FieldDefaults: get_default
import FieldMetadata: @units, units, @prior, prior, @description, description
import Flatten: flattenable

# Package to load and grid WOA data
using WorldOceanAtlasTools
