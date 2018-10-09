# Plots the trace of the optimizer in a 2d slice of the 3d parameter space,
# on top of a 2d contour of costs

using PyPlot, PyCall

function convergence_plot(results, q!)
    fig = figure(figsize=(8,6))
    ax = PyPlot.axes(yscale="log")
    out = Dict()
    
    # Determine number of plots to shift them
    n_lines = sum([haskey((results[k].trace[1]).metadata, "x") for k in keys(results)])
    counter_line = 0
    for k in mysort(keys(results))
        !haskey((results[k].trace[1]).metadata, "x") && continue
        λs = Optim.x_trace(results[k])
        y = [q!(λ) for λ in λs]
        ny = length(y)
        x = collect(1:ny) .+ counter_line / n_lines
        col = mycolor(k)
        #out[k] = semilogy(qvals, color = col, markevery=n_trace, marker="o", markersize=5, label=k)
        out[k] = step(x, y, color = col, markevery=[ny-1], marker="o", markersize=5, label=k, where="post")
        counter_line += 1
    end
    legend()
    #yaxis("log")
    xlim((1, 30))
    display(out)
    #savefig(out, "fig/convergence_plot_4box_model.pdf")
    return out
end

using Match
function mycolor(key)
    cmap = ColorMap("tab20c")
    println(key)
    println(cmapindex(key))
    println(cmap(cmapindex(key)))
    return cmap(cmapindex(key))
end

cmapindex(k) = @match k begin
    "Newton"              => 0
    "InteriorPointNewton" => 1
    "NewtonTrustRegion"   => 2
    "ConjugateGradient"   => 4
    "GradientDescent"     => 5
    "LBFGS"               => 6
    "SimulatedAnnealing"  => 8
    "NelderMead"          => 9
    "ParticleSwarm"       => 10
    _                     => 19
end

function mysort(keys)
    k_colors = [cmapindex(key) for key in keys]
    println(k_colors)
    return collect(keys)[sortperm(k_colors)]
end

convergence_plot(results, q!)