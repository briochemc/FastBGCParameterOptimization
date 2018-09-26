# Plots the trace of the optimizer in a 2d slice of the 3d parameter space,
# on top of a 2d contour of costs

using PyPlot, PyCall

function plot_figure_2(results, q!)

    clf()

    # 1) Determine the bounds of the traces
    λ_traces = vcat([(haskey((results[k].trace[1]).metadata, "x") ? Optim.x_trace(results[k]) : [λ₀]) for k in keys(results)]...)
    λ_extremas = extrema(hcat(λ_traces...), dims=2)

   #  x = range(λ_extremas[1][1], stop = λ_extremas[1][2], length = 100)
   #  y = range(λ_extremas[2][1], stop = λ_extremas[2][2], length = 80)
    x = range(-10.0, stop = 10.0, length = 100)
    y = range(-10.0, stop = 10.0, length = 80)
    z = [q!([x, y]) for y in y, x in x]

    filled_contour = contourf(x, y, z, levels = 0:0.001:0.03, extend = :max)
    cbar = colorbar(filled_contour)
    black_contours = contour(x, y, z, levels = 0:0.005:0.03, extend = :max, colors="black", linewidths=.5)
    cbar["add_lines"](black_contours)

    for k in keys(results)
        !haskey((results[k].trace[1]).metadata, "x") && continue
        xpath = hcat(Optim.x_trace(results[k])...)[1,:]
        ypath = hcat(Optim.x_trace(results[k])...)[2,:]
        col = rand(3)
        plt = plot(xpath, ypath, color = col, linestyle = "-", marker = "o", markersize = 5, markeredgecolor = "black", markerfacecolor = col, label = k)
    end

    legend()

    display(plt)

    savefig("Fig2.pdf")

end

plot_figure_2(results, q!)

results["Newton"].minimizer