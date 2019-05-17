# Load plotting package
using PyPlot

# Load data
using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
str_out = ""
jld_file = joinpath(path_to_package_root, "data", "TimerOutputs_data" * str_out * "" * ".jld2")
@load jld_file timers

translate_for_legend = Dict("f"=>"objective", "∇f"=>"gradient", "∇²f"=>"Hessian")

list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_∇ᵏf = ["objective", "gradient", "Hessian"]
list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_m = ["AF-1", "F-1", "DUAL", "CSD", "FD1", "HYPER", "FD2"]

tf = [timers[string(m)][string(m, "_f̂")]["time"] * 1e-9 for m in list_methods]
tg = [timers[string(m)][string(m, "_∇f̂")]["time"] * 1e-9 for m in list_methods]
th = [timers[string(m)][string(m, "_∇²f̂")]["time"] * 1e-9 for m in list_methods]

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
lighter_gray   = 0.85 * [1, 1, 1]
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
fig = figure("pyplot_barplot", figsize=(7,3))
subplots_adjust(left  = 0.12)  # the left side of the subplots of the figure
subplots_adjust(right = 0.99)
subplots_adjust(bottom = 0.15)
subplots_adjust(top = 0.98)
subplot(111, frameon=false)
ia = [7, 6, 4, 5, 3, 2] # indices of (a) plot
xa = 1:length(ia) # centers of bars
af = barh(xa, tf[ia], color = colf,
          align = "center",
          tick_label = list_m[ia])
ag = barh(xa, tg[ia], color = colg, align = "center", left = tf[ia])
ah = barh(xa, th[ia], color = colh, align = "center", left = tf[ia] + tg[ia])
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
legend((ag, ah), ("gradient","Hessian"), handlelength=1, handleheight=1)

eps_file = joinpath(path_to_package_root, "fig", "TimerOutputs_time.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "TimerOutputs_time.png")
savefig(png_file)


#=
subplot(211, frameon=false)
ia = [4, 5, 3, 2] # indices of (a) plot
xa = 1:length(ia) # centers of bars
af = barh(xa, tf[ia], color = colf,
          align = "center",
          tick_label = list_m[ia])
ag = barh(xa, tg[ia], color = colg, align = "center", left = tf[ia])
ah = barh(xa, th[ia], color = colh, align = "center", left = tf[ia] + tg[ia])
axis("tight")
grid("on", axis="x")
xlabel("")
ylabel("Method")
legend((af, ag, ah), ("obj","grad","hess"), handlelength=1, handleheight=1)

subplot(212, frameon=false)
ib = [7, 6, 2] # indices of (a) plot
xb = 1:length(ib) # centers of bars
bf = barh(xb, tf[ib], color = colf,
          align = "center",
          tick_label = list_m[ib])
bg = barh(xb, tg[ib], color = colg, align = "center", left = tf[ib])
bh = barh(xb, th[ib], color = colh, align = "center", left = tf[ib] + tg[ib])
axis("tight")
#PyPlot.title("(b)", loc="left")
grid("on", axis="x")
xlabel("Computation time (seconds)")
ylabel("Method")

PyPlot.suptitle("Partition of computation time")

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

df1 = timer_to_DataFrame(timers, [2,3,5,4])
df2 = timer_to_DataFrame(timers, [2,6,7])

orderf = [string(f) for f in list_∇ᵏf]
orderf = list_∇ᵏf
order1 = [string(m) for m in list_m[[2,3,5,4]]]
order2 = [string(m) for m in list_m[[2,6,7]]]

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


=#
