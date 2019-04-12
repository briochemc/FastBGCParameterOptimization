using JLD2, BenchmarkTools, VegaLite, DataFrames, Statistics

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results list_functions

list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_∇ᵏf = [:f, :∇f, :∇²f]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_m = ["F-1", "DUAL", "CSD", "FD1", "HYPER", "FD2"]

strings = string.(list_functions)
f̂_list = strings[[occursin("_f̂", s) for s in strings]]
∇f̂_list = strings[[occursin("_∇f̂", s) for s in strings]]
∇²f̂_list = strings[[occursin("_∇²f̂", s) for s in strings]]
f̂_to_f(x::String) = replace(replace(x, "f̂" => "f"), "F1" => "F-1")
f̂_to_f(x::Array) = [f̂_to_f(e) for e in x]


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
                :method => f̂_to_f(string(m, " ", fun)),
                :fgh => f̂_to_f(string(fun)),
                :time => results[f].times[] * 1e-9,
                :allocs => results[f].allocs[] * 2^-27
            )
        )
    end
    return df
end


df1 = results_to_df(results, f̂_list)
mymean(x) = mean(x)
function mymean(x::Array{String})
    if x[1] == "f"
        return "f"
    elseif x[1] == "∇f"
        return "∇f"
    elseif occursin("∇f", x[1])
        return "Analytical ∇f"
    else
        return "Analytical f"
    end
end
push!(df1, colwise(mymean, df1))
filter!(row -> row[:method] == "Analytical f", df1)

df2 = results_to_df(results, ∇f̂_list)
df2a = filter(row -> row[:method] ∈ ["DUAL ∇f", "CSD ∇f", "FD1 ∇f"], df2)
df2b = filter(row -> row[:method] ∉ ["DUAL ∇f", "CSD ∇f", "FD1 ∇f"], df2)
push!(df2a, colwise(mymean, df2a))
filter!(row -> row[:method] == "Analytical ∇f", df2a)
df2 = vcat(df2a, df2b)

df3 = results_to_df(results, ∇²f̂_list)

df = vcat(df1, df2, df3)
myround(x) = x ≥ 100 ? round(x) : round(x*10) / 10
df[:time] .= myround.(df[:time])

p = df |> @vlplot(
    #title={text="(a)", anchor="start"},
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

display(p)

# print to PDF
pdf_file = joinpath(path_to_package_root, "fig",  "BenchmarkTools_data" * str_out * ".pdf")
save(pdf_file, p)


