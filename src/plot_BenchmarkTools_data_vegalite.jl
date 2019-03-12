using JLD2, BenchmarkTools, VegaLite, DataFrames

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results


# Reorder keys to plot in my order
#     key       fgh       method
mykeys = [
    ("q!"    , "f"  , "f"         ),
    ("Dq!"   , "∇f" , "∇f"        ),
    ("ADq!"  , "∇f" , "Dual"      ),
    ("D2q!"  , "∇²f", "FLASH"     ),
    ("ADDq!" , "∇²f", "DUAL"      ),
    ("CSDDq!", "∇²f", "COMPLEX"   ),
    ("FDDq!" , "∇²f", "FINITEDIFF"),
    ("AD2q!" , "∇²f", "HYPERDUAL" )
]

#sorting = [
#    "f"
#    "∇f"
#    "Dual"
#    "FLASH"
#    "DUAL"
#    "COMPLEX"
#    "FINITEDIFF"
#    "HYPERDUAL"
#]

function results_to_df(results, mykeys)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0)
    )
    for (k, f, m) in mykeys
        df = push!(
            df,
            Dict(
                :method => m,
                :fgh => f,
                :time => round(results[k].times[] * 1e-9),
                :allocs => results[k].allocs[] * 2^-27
            )
        )
    end
    return df
end

df = results_to_df(results, mykeys)


p = df |> @vlplot(
    encoding={
        x={:time, title="Computation time (seconds)"},
        y={:method, sort=true},
        color={:fgh, title="", scale={scheme="set2"}},
        text={:time},
    },
    resolve={
        scale={
            y ="shared"
        }
    }
) +
@vlplot(:bar) + 
@vlplot(
    mark={
        :text,
        align=:left,
        baseline=:middle,
        dx=5
    }
)

display(p)

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

