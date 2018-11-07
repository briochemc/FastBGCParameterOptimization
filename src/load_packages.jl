# Main script to plot paper figures

# 1. Always load the packages
using Pkg
Pkg.activate(".")
using TransportMatrixTools
using SparseArrays, SuiteSparse, LinearAlgebra, Printf, DualNumbers
using JLD2, Calculus, Optim
# Packages for unitful parameters :)
using FieldDefaults, Flatten, FieldMetadata, Unitful
using Distributions
import FieldDefaults: get_default
import FieldMetadata: @units, units, @prior, prior, @description, description # 
import Flatten: flattenable


