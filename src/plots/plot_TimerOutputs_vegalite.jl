# Load plotting package
using VegaLite, DataFrames # for stacked bars

# Load data
using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "TimerOutputs_data" * str_out * "_katana" * ".jld2")
@load jld_file timers

translate_for_legend = Dict("f"=>"objective", "∇f"=>"gradient", "∇²f"=>"Hessian")

list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_∇ᵏf = [:f, :∇f, :∇²f]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_m = ["F-1", "DUAL", "CSD", "FD1", "HYPER", "FD2"]




function timer_to_DataFrame(timers::Dict, I)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0),
        ncalls = Array{Int64}(undef, 0)
    )
    for (m, m2) in zip(list_methods[I], list_m[I])
        t = timers[string(m)]
        for (∇ᵏf̂, ∇ᵏf) in zip(list_∇ᵏf̂, list_∇ᵏf)
            push!(
                df,
                Dict(
                    :method => string(m2),
                    :fgh => string(∇ᵏf),
                    :time => t[string(m, "_", ∇ᵏf̂)]["time"] * 1e-9,
                    :allocs => t[string(m, "_", ∇ᵏf̂)]["allocs"] * 2^-30,
                    :ncalls => t[string(m, "_", ∇ᵏf̂)]["ncalls"]
                )
            )
        end
    end
    return df
end

df1 = timer_to_DataFrame(timers, 1:4)
df2 = timer_to_DataFrame(timers, [1,5,6])

orderf = [string(f) for f in list_∇ᵏf]
order1 = [string(m) for m in list_m[1:4]]
order2 = [string(m) for m in list_m[[1,5,6]]]

ptime1 = df1 |> @vlplot(
    width=300,
    height=100,
    title={text="(a)", anchor="start"},
    mark={
        :bar,
    },
    x={"sum(time)", title="Computation time (seconds)", scale={domain=[0,3000]}},
    y={:method, sort=order1},
    color={:fgh, title="", scale={scheme="set2"}, legend={orient="top-right", sort=orderf}},
    order={field=:method, sort=orderf}
)

ptime2 = df2 |> @vlplot(
    width=300,
    height=75,
    title={text="(b)", anchor="start"},
    mark={
        :bar,
    },
    x={"sum(time)", title="Computation time (seconds)"},
    y={:method, sort=order2},
    color={:fgh, title="", scale={scheme="set2"}, legend={orient="top-right", sort=orderf}},
    order={field=:method, sort="count"}
)

ptime = [ptime1; ptime2]

display(ptime)

# print to PDF
pdf_file_time = joinpath(path_to_package_root, "fig", "TimerOutputs_time.pdf")
save(pdf_file_time, ptime)



