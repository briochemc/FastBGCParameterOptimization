using JLD2, Plots

# Load the benchmark data to plot
# TODO save the files with specific names to where it was run?
@load "data/methods_benchmark_OCIM1_imac.jld2"
# @load "data/methods_benchmark_data.jld2"
msbd = methods_benchmark_data

# TODO save the list_methods with the convergence results
# TODO save the final params and cost to compare them too
list_methods = [
    ( "D2q", :q!,  :Dq!,   :D2q!)
    ("AD2q", :q!, :ADq!,  :AD2q!)
    ("CSDq", :q!,  :Dq!, :CSDDq!)
    ("FDDq", :q!,  :Dq!,  :FDDq!)
]

fields_with_time = [
    "q",
    "Dq",
    "FDq",
    "CSDq",
    "ADq",
    "D2q",
    "FD2q",
    "FDDq",
    "CSDDq",
    "ADDq",
    "AD2q",
    "factorize",
    "backsolve"
]

start_times = Dict()
for method in list_methods
    #for run_str in ["_1strun", "_2ndrun"]
    for run_str in ["_1strun"]
        field_name = method[1] * run_str
        mbd = msbd[field_name]
        #println(field_name)
        #println(minimum([first(mbd[f][1]) for f in fields_with_time if has_timer(mbd, f)]))
        #println([msbd[field_name][ftime][1] for ftime in fields_with_time if ~isempty(msbd[field_name][ftime][1])])
        push!(start_times, field_name => minimum([first(mbd[f][1]) for f in fields_with_time if has_timer(mbd, f)]))
    end
end

has_timer(mbd, f) = haskey(mbd, f) && ~isempty(mbd[f][1])


end_costs = Vector{Float64}()
for method in list_methods
    for run_str in ["_1strun"]
    #for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        push!(end_costs, minimum(methods_benchmark_data[field_name]["qvalues"]))
    end
end
using Statistics
end_cost = mean(end_costs)

p = plot()
for (i, method) in enumerate(list_methods)
    for run_str in ["_1strun"]
        field_name = method[1] * run_str
        x = (methods_benchmark_data[field_name]["q"][2] .- start_times[field_name]) ./ 60
        y = abs.(methods_benchmark_data[field_name]["qvalues"] .- end_cost)
        plot!(p, x, y, yaxis = :log, label = field_name)
        xlabel!(p, "time (min)")
        ylabel!(p, "cost function error")
        display(p)
    end
end
display(p)


###
#
# Start Checking time spent in each phase
#
###
markers = [:o, :d, :s, :v]

p = plot()
for (i, method) in enumerate(list_methods)
    field_name = method[1] * "_1strun"
    mbd = msbd[field_name]
    for f in fields_with_time
        f == "q" || isempty(mbd[f][1]) ? continue : nothing
        println(f)
        tictocs = mbd[f]
        label = method[1] * " - " * f
        plot!(p, cumsum(tictocs[2] - tictocs[1])/60, axis = :log, label = label, marker = markers[i])
        display(p)
    end
end


####
#
# End Checking time spent in each phase
#
####

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
    x, y, oldv = Vector{Float64}(), Vector{Float64}(), NaN
    for (tb, ta, v) in zip(tbefore, tafter, vals)
        push!(x, tb)
        push!(y, oldv)
        push!(x, ta)
        push!(y, v)
        oldv = v
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

