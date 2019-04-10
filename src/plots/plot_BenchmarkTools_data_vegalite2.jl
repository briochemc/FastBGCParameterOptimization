using JLD2, BenchmarkTools, VegaLite, DataFrames

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results list_functions

strings = string.(list_functions)
f̂_list = strings[[occursin("_f̂", s) for s in strings]]
∇f̂_list = strings[[occursin("_∇f̂", s) for s in strings]]
∇²f̂_list = strings[[occursin("_∇²f̂", s) for s in strings]]

function results_to_df(results, list)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0)
    )
    for f in list
        m, fun = split(f, "_")
        df = push!(
            df,
            Dict(
                :method => m,
                :fgh => fun,
                :time => myround(results[f].times[] * 1e-9),
                :allocs => results[f].allocs[] * 2^-27
            )
        )
    end
    return df
end

myround(x) = x ≥ 100 ? round(x) : round(x*10) / 10


df1 = results_to_df(results, f̂_list)
df2 = results_to_df(results, ∇f̂_list)
df3 = results_to_df(results, ∇²f̂_list)

p1 = df1 |> @vlplot(
    title={text="(a)", anchor="start"},
    encoding={
        x={:time, title="Computation time (seconds)"},
        y={:method, sort=true},
        color={:fgh, title="", scale={scheme="set2"}, legend={orient="top-right"}},
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

p2 = df2 |> @vlplot(
    title={text="(b)", anchor="start"},
    encoding={
        x={:time, title="Computation time (seconds)", scale={domain=[0,2500]}},
        y={:method, sort=true},
        color={:fgh, title="", scale={scheme="set2"}, legend={orient="top-right"}},
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

p3 = df3 |> @vlplot(
    title={text="(b)", anchor="start"},
    encoding={
        x={:time, title="Computation time (seconds)", scale={domain=[0,2500]}},
        y={:method, sort=true},
        color={:fgh, title="", scale={scheme="set2"}, legend={orient="top-right"}},
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

p = [p1; p2; p3]

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

