using JLD2, Plots

# Load the Cassette profile data to plot
# TODO save the files with specific names to where it was run?
# @load "data/Cassette_profile_OCIM1_imac.jld2"
@load "data/Cassette_profile_OCIM1_katana.jld2"
# @load "data/Cassette_profile.jld2"
K7profiles = Cassette_profile

# TODO save the list_methods with the convergence results
# TODO save the final params and cost to compare them too
list_methods = [
    ( "D2q", :q!,  :Dq!,   :D2q!)
    ("AD2q", :q!, :ADq!,  :AD2q!)
    ("CSDq", :q!,  :Dq!, :CSDDq!)
    ("FDDq", :q!,  :Dq!,  :FDDq!)
]

timed_functions = [
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
has_timer(K7profile, f) = haskey(K7profile, f) && ~isempty(K7profile[f][1])
for method in list_methods
    for run_str in ["_1strun", "_2ndrun"]
    #for run_str in ["_1strun"]
        field_name = method[1] * run_str
        K7profile = K7profiles[field_name]
        #println(field_name)
        #println(minimum([first(K7profile[f][1]) for f in timed_functions if has_timer(K7profile, f)]))
        #println([K7profiles[field_name][ftime][1] for ftime in timed_functions if ~isempty(K7profiles[field_name][ftime][1])])
        push!(start_times, field_name => minimum([first(K7profile[f][1]) for f in timed_functions if has_timer(K7profile, f)]))
    end
end



end_costs = Vector{Float64}()
for method in list_methods
    #for run_str in ["_1strun"]
    for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        push!(end_costs, minimum(Cassette_profile[field_name]["qvalues"]))
    end
end
using Statistics
end_cost = mean(end_costs)

p = plot()
for (i, method) in enumerate(list_methods)
    for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        x = (Cassette_profile[field_name]["q"][2] .- start_times[field_name]) ./ 60
        y = abs.(Cassette_profile[field_name]["qvalues"] .- end_cost)
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
    field_name = method[1] * "_2ndrun"
    K7profile = K7profiles[field_name]
    println(method[1], ":")
    for f in timed_functions
        isempty(K7profile[f][1]) ? continue : nothing
        tictocs = K7profile[f]
        label = method[1] * " - " * f
        #plot!(p, cumsum(tictocs[2] - tictocs[1])/60, axis = :log, label = label, marker = markers[i])
        plot!(p, cumsum(tictocs[2] - tictocs[1])/60, label = label, marker = markers[i])
        display(p)
        println("- ", f, ": ", round(sum(tictocs[2] - tictocs[1])/60), " min")
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
        push!(min_ys, minimum(Cassette_profile[field_name]["qvalues"]))
    end
end
best_y = minimum(min_ys)

# plot vs timer
p = plot() ;
for (i, method) in enumerate(list_methods)
    for run_str in ["_1strun", "_2ndrun"]
        field_name = method[1] * run_str
        x = Cassette_profile[field_name]["qtimes_before"] / 60
        y = Cassette_profile[field_name]["costs"]
        y = y .- best_y
        y .= map(y -> y ≤ 0 ? NaN : y, y)
        plot!(p, x, y, yaxis = :log, label = field_name)
        xlabel!(p, "time (min)")
        ylabel!(p, "error in cost function")
        display(p)
    end
end
display(p)
# savefig("fig/Cassette_profile_logDeltaq_vs_time.pdf")

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
# savefig("fig/Cassette_profile_logDeltaq_vs_nfact.pdf")


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

