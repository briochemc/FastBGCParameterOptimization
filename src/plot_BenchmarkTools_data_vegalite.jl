using JLD2, BenchmarkTools, VegaLite

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results

for k in keys(results)
    println(results[k])
end


# Reorder keys to plot in my order
#     key       fgh       method
mykeys = [
    ("q!"    , "q"  , "q"     ),
    ("Dq!"   , "dq" , "dq"     ),
    ("D2q!"  , "d²q", "D2q"  ),
    ("FDDq!" , "d²q", "FDDq" ),
    ("ADDq!" , "d²q", "ADDq" ),
    ("CSDDq!", "d²q", "CSDDq")
]


function results_to_df(results, mykeys)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0)
    )
    for (k, f, m) in mykeys
        println(typeof(results[k].times[]))
        println(results[k].times)
        df = push!(
            df,
            Dict(
                :method => m,
                :fgh => f,
                :time => results[k].times[] * 1e-9,
                :allocs => results[k].allocs[] * 2^-27
            )
        )
    end
    return df
end

df = results_to_df(results, mykeys)


p = df |> @vlplot(
    :bar,
    x={:time, title="Computation time (seconds)"},
    y={:method, sort=["q", "dq", "D2q", "FDDq", "ADDq", "CSDDq"]},
    color={:fgh, title="", scale={scheme="set2"}, sort=["q", "dq", "d²q"]},
    order={field=:method, sort=["q", "dq", "d²q"]}
)

# print to PDF
pdf_file = joinpath(path_to_package_root, "fig",  "BenchmarkTools_data" * str_out * ".pdf")
save(pdf_file, p)





#p = plot()
#for (i, k) in enumerate(mykeys)
#    println(float(i), " - ", results[k].times[] * 1e-9)
#    x = results[k].times[] * 1e-9 # time in seconds
#    mycolor = mycolors[2]
#    k == "q!" ? mycolor = mycolors[3] : nothing
#    k == "Dq!" ? mycolor = mycolors[1] : nothing
#    plot!(p, (i, x), seriestype=:bar, label=k[1:end-1], color=mycolor)
#end
#plot!(p, yticks=collect(0:60:500))
#plot!(p, ylabel="Time (seconds)")
#plot!(p, xticks=(1:6,map(x->x[1:end-1], mykeys)))
#plot!(p, legend=nothing)
#display(p)
#
#savefig(p, "fig/BenchmarkTools_data_katana2.pdf")

