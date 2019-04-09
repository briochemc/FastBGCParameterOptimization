using JLD2, BenchmarkTools

#@load "data/BenchmarkTools_data.jld2"
@load "data/BenchmarkTools_data_katana.jld2"

for k in keys(results)
    println(results[k])
end

using Plots
gr()

using ColorBrewer # for better colors

# Reorder keys to plot in my order
mykeys = [
    "q!"
    "Dq!"
    "D2q!"
    "FDDq!"
    "ADDq!"
    "CSDDq!"
]

mycolors = ColorBrewer.palette("Set2", 6)

p = plot()
for (i, k) in enumerate(mykeys)
    println(float(i), " - ", results[k].times[] * 1e-9)
    x = results[k].times[] * 1e-9 # time in seconds
    mycolor = mycolors[2]
    k == "q!" ? mycolor = mycolors[3] : nothing
    k == "Dq!" ? mycolor = mycolors[1] : nothing
    plot!(p, (i, x), seriestype=:bar, label=k[1:end-1], color=mycolor)
end
plot!(p, yticks=collect(0:60:500))
plot!(p, ylabel="Time (seconds)")
plot!(p, xticks=(1:6,map(x->x[1:end-1], mykeys)))
plot!(p, legend=nothing)
display(p)

savefig(p, "fig/BenchmarkTools_data_katana2.pdf")

