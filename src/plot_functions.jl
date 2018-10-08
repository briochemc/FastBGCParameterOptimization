# Plots the trace of the optimizer in a 2d slice of the 3d parameter space,
# on top of a 2d contour of costs

using PyPlot, PyCall

function convergence_plot(results, q!)
    clf()
    out = nothing
    for k in keys(results)
        !haskey((results[k].trace[1]).metadata, "x") && continue
        x_trace = Optim.x_trace(results[k])
        qvals = [q!(λ) for λ in x_trace]
        col = mycolor(k)
        out = semilogy(qvals, color = col)
        semilogy(length(qvals), qvals[end], color = col, marker="o", markersize=5)
    end
    legend(keys(results))
    xlim((1, 10))
    display(out)
    savefig(out, "fig/convergence_plot_4box_model.pdf")
    return out
end

using Match
function mycolor(key)
    cmap = ColorMap("tab20c")
    @match key begin
        "Newton"              => cmap(0)
        "InteriorPointNewton" => cmap(1)
        "NewtonTrustRegion"   => cmap(2)
        "ConjugateGradient"   => cmap(4)
        "GradientDescent"     => cmap(5)
        "LBFGS"               => cmap(6)
        "SimulatedAnnealing"  => cmap(8)
        "NelderMead"          => cmap(9)
        "ParticleSwarm"       => cmap(10)
        _                     => cmap(19)
    end
end

convergence_plot(results, q!)
