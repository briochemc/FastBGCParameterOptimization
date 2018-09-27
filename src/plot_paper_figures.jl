# Main script to plot paper figures

# 1. Always load the packages
using Pkg
Pkg.activate(".")
using TransportMatrixTools
using SparseArrays, SuiteSparse, LinearAlgebra, Printf, DualNumbers
using JLD2, Parameters, Optim
using Flatten
import Flatten: flattenable



# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)
# include("setup_OCIM1.1.jl")
include("setup_4BoxModel.jl")
# The above must define everything needed by the functions

# load biogeochmistry functions
include("bgc_functions.jl")

# load cost functions
include("cost_functions.jl")


# # print("\nOptimizing (Q Newton)...")
# # @benchmark print_results(optimize(Qwrap, λ₀, Newton()))
# # print("\nOptimizing (Newton with FD gradient and Hessian)...")
# # @benchmark print_results(optimize(slowQ, FDQ!, FD2Q!, λ₀, Newton()))
# # print("\nOptimizing (Newton with gradient and Hessian)...")

function print_results(results)
    λopt = results.minimizer
    popt = λ2p(λopt)
    paraprint(popt)
    qopt = q!(popt)
    print_cost(qopt)
    printstyled("  gives cost q = $(@sprintf("%.2g",qopt))\n", color=:red)
    @show results
  #println("    which means ≈ $(@sprintf("%.2g",(100*sqrt(qopt))))\% RMS")
  #println("      And steady state a∞ = $(newton_solver(popt))")
  #println(results)
end

# Do those once for the benchmark
using BenchmarkTools
q!(λ₀)
Dq!(λ₀)
D2q!(λ₀)

opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = true, x_tol = 1e-3)
results = Dict()
# Hessian Required
  # 1) Newton
  λshift = 0.02 * ones(npopt)
  λ₀ = λ₀ + λshift
@btime results["Newton"] = optimize(q!, Dq!, D2q!, λ₀, Newton(), opt)
# print_results(results["Newton"])
  # 2) Interior Point Newton
  λ₀ = λ₀ + λshift
df = TwiceDifferentiable(q!, Dq!, D2q!, λ₀)
lx = Float64[]; ux = Float64[]
dfc = TwiceDifferentiableConstraints(lx, ux)
@btime results["InteriorPointNewton"] = optimize(df, dfc, λ₀, IPNewton(), opt)
# print_results(results["InteriorPointNewton"])
  # 3) Newton Trust Region
  λ₀ = λ₀ + λshift
@btime results["NewtonTrustRegion"] = optimize(q!, Dq!, D2q!, λ₀, NewtonTrustRegion(), opt)
# print_results(results["NewtonTrustRegion"])
# Gradient Required
  # 4) Conjugate Gradient - avoid singular expression by normalizing by gradient at starting point
  λ₀ = λ₀ + λshift
@btime results["ConjugateGradient"] = optimize(q!, Dq!, λ₀, ConjugateGradient(; P = norm(Dq!(λ₀)) * Matrix(I, npopt, npopt)), opt)
# print_results(results["ConjugateGradient"])
  # 5) Gradient Descen - avoid singular expression by normalizing by gradient at starting point
  λ₀ = λ₀ + λshift
@btime results["GradientDescent"] = optimize(q!, Dq!, λ₀, GradientDescent(; P = norm(Dq!(λ₀)) * Matrix(I, npopt, npopt)), opt)
# print_results(results["GradientDescent"])
  # 6) LBFGS - avoid singular expression by normalizing by gradient at starting point
  λ₀ = λ₀ + λshift
@btime results["LBFGS"] = optimize(q!, Dq!, λ₀, LBFGS(; P = norm(Dq!(λ₀)) * Matrix(I, npopt, npopt)), opt)
# print_results(results["LBFGS"])
# Gradient Free
  # 7) Particle Swarm - I struggle with this one...
# lλ = [10.0, -11.0]
# uλ = [16.0, -7.0]
# results["ParticleSwarm"] = optimize(q!, λ₀, ParticleSwarm(; lower = lλ, upper = uλ), opt)
# print_results(results["ParticleSwarm"])
  # 8) Simulated Annealing
  λ₀ = λ₀ + λshift
@btime results["SimulatedAnnealing"] = optimize(q!, λ₀, SimulatedAnnealing(), opt)
# print_results(results["SimulatedAnnealing"])
  # 9) NelderMead
  λ₀ = λ₀ + λshift
@btime results["NelderMead"] = optimize(q!, λ₀, NelderMead(), opt)
print_results(results["NelderMead"])


# include("plot_figure_1.jl")
#include("plot_figure_2.jl")


# λopt = results.minimizer
# popt = λ2p(λopt)

# xsol = solvef0NewtonChordShamanskii!(x -> f(x, popt), x -> fJac(x, popt), x₀, nrm, τstop, true)
