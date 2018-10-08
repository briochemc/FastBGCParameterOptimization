# Main script to plot paper figures

# 1. Always load the packages
include("load_packages.jl")

# Setup OCIM1.1 or toy model (comment/uncomment to use the one you need for now)
# include("setup_OCIM1.1.jl")
include("setup_4BoxModel.jl")
# The above must define everything needed by the functions

# load biogeochmistry parameters
include("bgc_parameters.jl")

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
    show(popt)
    qopt = q!(popt)
    print_cost(qopt)
    printstyled("  gives cost q = $(@sprintf("%.2g",qopt))\n", color=:red)
    @show results
  #println("    which means ≈ $(@sprintf("%.2g",(100*sqrt(qopt))))\% RMS")
  #println("      And steady state a∞ = $(newton_solver(popt))")
  #println(results)
end

using JLD2
if false # write as true if you want to run the optimizers
    include("run_optimizers.jl")
    @save "fig/4BoxModel_results.jld2" results
else
    @load "fig/4BoxModel_results.jld2" results
end

# include("plot_figure_1.jl")
# include("plot_figure_2.jl")
# include("plot_figure_3.jl")
include("plot_functions.jl")
