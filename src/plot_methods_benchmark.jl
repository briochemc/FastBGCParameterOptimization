using JLD2, Plots

# Load the benchmark data to plot
# TODO save the files with specific names to where it was run?
@load "data/methods_benchmark_OCIM1_imac.jld2"

# TODO save the list_methods with the convergence results
# TODO save the final params and cost to compare them too
list_methods = [
    ( "D2q", :q!,  :Dq!,   :D2q!)
    ("AD2q", :q!, :ADq!,  :AD2q!)
    ("CSDq", :q!,  :Dq!, :CSDDq!)
    ("FDDq", :q!,  :Dq!,  :FDDq!)
]

min_ys = Vector{Float64}()
for (i, method) in enumerate(list_methods)
    for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        push!(min_ys, minimum(methods_benchmark_data[field_name]["qvalues"]))
    end
end
best_y = minimum(min_ys)

# plot vs timer
p = plot() ;
for (i, method) in enumerate(list_methods)
    for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        x = methods_benchmark_data[field_name]["qtimes_before"] / 60
        y = methods_benchmark_data[field_name]["costs"]
        y = y .- best_y
        y .= map(y -> y ≤ 0 ? NaN : y, y)
        plot!(p, x, y, yaxis = :log, label = field_name)
        xlabel!(p, "time (min)")
        ylabel!(p, "error in cost function")
        display(p)
    end
end
display(p)
# savefig("fig/methods_benchmark_logDeltaq_vs_time.pdf")

function mystep(tbefore, tafter, vals)
    x, y = Vector{Float64}(), Vector{Float64}()
    for (tb, ta, v) in zip(tbefore, tafter, vals)
        push!(x, tb)
        push!(y, v)
        push!(x, ta)
        push!(y, v)
    end
    return x, y
end


# plot vs factorizations
p = plot() ;
for (i, method) in enumerate(list_methods)
    run_str = "_2ndrun"
    field_name = method[1] * run_str
    x = convergence_results[field_name]["facts"]
    y = convergence_results[field_name]["costs"]
    y = y .- best_y
    y .= map(y -> y ≤ 0 ? NaN : y, y)
    plot!(p, x, y, yaxis = :log, label = method[1])
    xlabel!(p, "number of calls to `factorize`")
    ylabel!(p, "error in cost function")
    display(p)
end
display(p)
# savefig("fig/methods_benchmark_logDeltaq_vs_nfact.pdf")


# plot vs factorizations
p = plot() ;
for (i, method) in enumerate(list_methods)
    run_str = "_2ndrun"
    field_name = method[1] * run_str
    x = convergence_results[field_name]["facts"]
    y = convergence_results[field_name]["times"]
    plot!(p, x, y)
    display(p)
end
display(p)

# # Now plot the results using `Plots.jl`
# using Plots
# 
# # Convergence vs time
# timer = tape.metadata.ftimer .- tape.metadata.ftimer[1]
# y = tape.metadata.fvalues
# p1 = plot(timer, y, yaxis = :log)
# xlabel!("computing time (ms)")
# ylabel!("f")
# legend!("")
# 
# # Convergence vs number of calls
# counter = tape.metadata.fcounter
# p2 = plot(counter, y, yaxis = :log)
# xlabel!("number of f calls")
# ylabel!("f")
# 
# # Combine subplotds into single figure
# plot(p1, p2, layout = (2, 1))
#

