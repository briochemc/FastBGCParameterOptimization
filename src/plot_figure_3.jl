# Plots the trace of the optimizer in a 2d slice of the 3d parameter space,
# on top of a 2d contour of costs

using PyPlot, PyCall

function plot_figure_3(results, q!)

    clf()

    # 1) Determine the bounds of the traces
    λ_traces = vcat([Optim.x_trace(results[k]) for k in keys(results)]...)
    λ_extremas = extrema(hcat(λ_traces...), dims=2)

   #  x = range(λ_extremas[1][1], stop = λ_extremas[1][2], length = 100)
   #  y = range(λ_extremas[2][1], stop = λ_extremas[2][2], length = 80)

    for k in keys(results)
        x_trace = Optim.x_trace(results[k])
        y = [q!(λ) for λ in x_trace]
        col = rand(3)
        plt = plot(y, color = col, linestyle = "-", marker = "o", markersize = 5, markeredgecolor = "black", markerfacecolor = col, label = k)
    end

    legend()

    display(plt)

    savefig(plt, "Fig3.pdf")

end

plot_figure_3(results, q!)

results["Newton"].minimizer