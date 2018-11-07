# Main script to plot paper figures

# 1. Always load the packages
using Pkg
Pkg.activate(".")
using TransportMatrixTools
using SparseArrays, SuiteSparse, LinearAlgebra, Printf, DualNumbers
using JLD2, Calculus, Optim
# Packages for unitful parameters :)
using FieldDefaults, Flatten, FieldMetadata, Unitful
import FieldDefaults: get_default
import FieldMetadata: @units, units # use the metadata of units already in it
import Flatten: flattenable


