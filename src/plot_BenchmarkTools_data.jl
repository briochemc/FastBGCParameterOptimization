using JLD2, BenchmarkTools

@load "data/BenchmarkTools_data.jld2"

for k in keys(results)
    println(results[k])
end

using Plots

for (i, k) in enumerate(keys(results))
    println(i, " - ", results[k].times[] * 1e-9 / 60)
    bar!(i, results[k].times[] * 1e-9 / 60, label=k)
end

bar(1,2)

