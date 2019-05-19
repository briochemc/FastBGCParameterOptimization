using PyPlot, JLD2, BenchmarkTools, DataFrames, StatsBase

#@load "data/BenchmarkTools_data.jld2"
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
str_out = ""
jld_file = joinpath(path_to_package_root, "data", "BenchmarkTools_data" * str_out * ".jld2")
@load jld_file results list_functions

list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_m = ["F-1", "DUAL", "CSD", "FD1", "HYPER", "FD2"]



#===========================================
Data into DataFrame 
===========================================#

strings = string.(list_functions)[4:end] # Remove AF1 for fairness in the paper
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
df2 = results_to_df(results, ∇f̂_list)
df3 = results_to_df(results, ∇²f̂_list)

myround(x) = x ≥ 100 ? round(x) : round(x*10) / 10
tf = myround.(df1[:time])
tg = myround.(df2[:time])
th = myround.(df3[:time])

#===========================================
Julia colors 
===========================================#

darker_blue    = [0.251, 0.388, 0.847]
lighter_blue   = [0.4  , 0.51 , 0.878]
darker_purple  = [0.584, 0.345, 0.698]
lighter_purple = [0.667, 0.475, 0.757]
darker_green   = [0.22 , 0.596, 0.149]
lighter_green  = [0.376, 0.678, 0.318]
darker_red     = [0.796, 0.235, 0.2  ]
lighter_red    = [0.835, 0.388, 0.361]
lighter_gray   = 0.90 * [1, 1, 1]
darker_gray    = 0.75 * [1, 1, 1]
using Colors
# function to transform RGB (from 0 to 1) color array into hex string
function convert_rgb_to_hex(aRGB)
    r, g, b = aRGB
    return "#" * hex(RGB(r, g, b))
end
colf = convert_rgb_to_hex(lighter_gray)
colg = convert_rgb_to_hex(darker_gray)
colh = convert_rgb_to_hex(darker_purple)

#===========================================
Plot 
===========================================#

fig = figure("pyplot_barplot", figsize=(5,6))
subplots_adjust(left  = 0.16)  # the left side of the subplots of the figure
subplots_adjust(right = 0.99)
subplots_adjust(bottom = 0.08)
subplots_adjust(top = 0.98)
subplot(111, frameon=false)
ia = 6:-1:1 # indices of plot
xa = 3 * (1:length(ia)) # centers of bars
w = 0.9 # bar height
af = barh(xa .+ w, tf[ia], w, color = colf, align = "center")
ag = barh(xa, tg[ia], w, color = colg, align = "center", tick_label = list_m[ia])
ah = barh(xa .- w, th[ia], w, color = colh, align = "center")
# Optional text bits
#fs = 5 # font size of text bits
#for (x, i) in enumerate(ia[1:end-1])
#    text(tf[i] + tg[i] + th[i], xa[x], string(Int(round(th[i]))), horizontalalignment="right", color="white", verticalalignment="center", fontsize=fs)
#    text(tf[i] + tg[i], xa[x], string(Int(round(tg[i]))), horizontalalignment="right", verticalalignment="center", fontsize=fs)
#end
#text(tf[2] + tg[2] + th[2], xa[end], string(Int(round(th[2]))), horizontalalignment="left", color=colh, verticalalignment="center", fontsize=fs)
#text(tf[2] + tg[2], xa[end], string(Int(round(tg[2]))), horizontalalignment="right", verticalalignment="center", fontsize=fs)
axis("tight")
grid("on", axis="x")
xlabel("Computation time (seconds)")
ylabel("Method")
xlim((0,3250))
legend((af, ag, ah), ("objective","gradient","Hessian"), handlelength=1, handleheight=1)

eps_file = joinpath(path_to_package_root, "fig", "BenchmarkTools_data.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "BenchmarkTools_data.png")
savefig(png_file)


#=
tf = [timers[string(m)][string(m, "_f̂")]["time"] * 1e-9 for m in list_methods]
tg = [timers[string(m)][string(m, "_∇f̂")]["time"] * 1e-9 for m in list_methods]
th = [timers[string(m)][string(m, "_∇²f̂")]["time"] * 1e-9 for m in list_methods]
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
=#
