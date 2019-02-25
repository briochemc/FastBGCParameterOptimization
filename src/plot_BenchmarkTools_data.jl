using JLD2, BenchmarkTools

#@load "data/BenchmarkTools_data.jld2"
@load "data/BenchmarkTools_data_katana.jld2"

for k in keys(results)
    println(results[k])
end

using Plots
gr()

markers = Dict(
    "D2q!" => :star5,
    "FDDq!" => :rect,
    "ADDq!" => :diamond,
    "CSDDq!" => :circle,
    "q!" => :+,
    "Dq!" => :x
)

# Reorder keys to plot in my order
mykeys = [
    "q!"
    "Dq!"
    "D2q!"
    "FDDq!"
    "ADDq!"
    "CSDDq!"
]

marker_colors = Dict(
    "D2q!" =>   myRGB( 31,120,180),
    "FDDq!" =>  myRGB(166,206,227),
    "ADDq!" =>  myRGB(178,223,138),
    "CSDDq!" => myRGB( 51,160, 44),
    "q!" =>  myRGB(0,0,0),
    "Dq!" => myRGB(0,0,0)
)

mylabels = Dict(
    "D2q!" =>   "This paper",
    "FDDq!" =>  "Finite diff",
    "ADDq!" =>  "Dual numbers",
    "CSDDq!" => "Complex step",
    "q!" =>  "q",
    "Dq!" => "Dq"
)

marker_sizes = Dict(
    "D2q!" =>   10.0,
    "FDDq!" =>  5.0,
    "ADDq!" =>  5.0,
    "CSDDq!" => 5.0,
    "q!" =>  5.0,
    "Dq!" => 5.0
)

myRGB(r,g,b) = RGB(r/255, g/255, b/255)


p = plot()
for (i, k) in enumerate(mykeys)
    println(float(i), " - ", results[k].times[] * 1e-9)
    x = results[k].times[] * 1e-9 # time in seconds
    y = results[k].memory * 1e-9 # memory in GB
    plot!(p, (x,y), line=nothing, label=mylabels[k], marker=markers[k], color=marker_colors[k], markersize=marker_sizes[k])
end
plot!(p, legend=:topleft, legendline=nothing)
plot!(p, xticks=collect(0:60:500))
plot!(p, xlabel="Time (seconds)")
plot!(p, ylabel="Allocations (GB)")
display(p)

savefig(p, "fig/BenchmarkTools_data_katana.pdf")

