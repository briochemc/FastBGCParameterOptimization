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
using DualNumbers, Calculus

# Loading and saving data
using JLD2, MAT

# Optimization package
using Optim

# Packages for parameters
using FieldDefaults, Flatten, FieldMetadata, Unitful, Distributions
import FieldDefaults: get_default
import FieldMetadata: @units, units, @prior, prior, @description, description
import Flatten: flattenable


