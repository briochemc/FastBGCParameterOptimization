using JLD2, BenchmarkTools, VegaLite, DataFrames

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results


# Reorder keys to plot in my order
#     key       fgh       method
mykeys = [
    ("f̂!"        , "f"  , "Analytical f" ),
    ("∇f̂!"       , "∇f" , "Analytical ∇f"), #┐
    ("OF1_∇f̂!"   , "∇f" , "F-1 ∇f"       ), #│
    ("F1_∇f̂!"    , "∇f" , "old F-1 ∇f"   ), #│
    ("HYPER_∇f̂!" , "∇f" , "HYPER ∇f"     ), #├─ gradients
    ("FD2_∇f̂!"   , "∇f" , "FD2 ∇f"       ), #┘
    ("OF1_∇²f̂!"  , "∇²f", "F-1 ∇²f"      ), #┐
    ("F1_∇²f̂!"   , "∇²f", "old F-1 ∇²f"  ), #│
    ("F0_∇²f̂!"   , "∇²f", "old F-0 ∇²f"  ), #│
    ("DUAL_∇²f̂!" , "∇²f", "DUAL ∇²f"     ), #│
    ("CSD_∇²f̂!"  , "∇²f", "CSD ∇²f"      ), #│
    ("FD1_∇²f̂!"  , "∇²f", "FD1 ∇²f"      ), #├─ Hessians
    ("HYPER_∇²f̂!", "∇²f", "HYPER ∇²f"    ), #│
    ("FD2_∇²f̂!"  , "∇²f", "FD2 ∇²f"      )  #┘
]

#sorting = [
#    "f"
#    "∇f"
#    "Dual"
#    "F-zero"
#    "DUAL"
#    "COMPLEX"
#    "FINITEDIFF"
#    "HYPERDUAL"
#]

function results_to_df(results, mykeys, m_list)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0)
    )
    for (k, f, m) in mykeys 
        if m ∈ m_list
            df = push!(
                df,
                Dict(
                    :method => m,
                    :fgh => f,
                    :time => myround(results[k].times[] * 1e-9),
                    :allocs => results[k].allocs[] * 2^-27
                )
            )
        end
    end
    return df
end

myround(x) = x ≥ 100 ? round(x) : round(x*10) / 10

m_list1 = [
    "Analytical f"
    "Analytical ∇f"
    "F-1 ∇f"
#    "old F-1 ∇f"
    "HYPER ∇f"
    "FD2 ∇f"
    "F-1 ∇²f"
#    "old F-1 ∇²f"
#    "old F-0 ∇²f"
    "DUAL ∇²f"
    "CSD ∇²f"
    "FD1 ∇²f"
]
m_list2 = [
    "F-1 ∇²f"
    "HYPER ∇²f"
    "FD2 ∇²f"
]
df1 = results_to_df(results, mykeys, m_list1)
df2 = results_to_df(results, mykeys, m_list2)

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

p = [p1; p2]

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

