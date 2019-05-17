using PyPlot #



#=
df = DataFrame(
    method = Array{String}(undef, 0),
    time = Array{Float64}(undef, 0),
    iteration = Array{Int64}(undef, 0),
    qval = Array{Float64}(undef, 0),
    normgradq = Array{Float64}(undef, 0)
)
function mypush!(df, method_name, time_iter_q_normgradq)
    start_time = time_iter_q_normgradq[1, 1]
    for i in 2:size(time_iter_q_normgradq, 1)
        time_iter_q_normgradq[i, 1] -= start_time
        time_iter_q_normgradq[:, 1] .= max.(time_iter_q_normgradq[:, 1], 100.0)
        push!(df, [method_name; time_iter_q_normgradq[i, :]])
    end
    return df
end
=#

list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
list_m = ["AF-1", "F-1", "DUAL", "CSD", "FD1", "HYPER", "FD2"]

x = Dict()
y = Dict()

println("The data being plotted was copy-pasted from Katana output!")

#===========================================
Data 
===========================================#
#    time                   iteration      q(λ)           |Dq(λ)|
m = "F-1"
time_iter_q_normgradq = [
    1.557751180373945e9        -1         NaN             NaN
    1.5577513210150471e+09     0     6.898228e-01     5.399008e-03
    1.5577513850427470e+09     1     6.826155e-01     5.289488e-03
    1.5577514285317619e+09     2     6.686786e-01     5.073548e-03
    1.5577514612624371e+09     3     6.422890e-01     4.696714e-03
    1.5577514939592800e+09     4     5.961951e-01     3.871112e-03
    1.5577515265045810e+09     5     5.304139e-01     2.181979e-03
    1.5577515590177331e+09     6     5.000000e-01     5.799814e-06
    1.5577515800671041e+09     7     4.999999e-01     8.587957e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "DUAL"
time_iter_q_normgradq = [
    1.557751580222345e9        -1         NaN             NaN
    1.5577518531085820e+09     0     6.898228e-01     5.399008e-03
    1.5577520433873570e+09     1     6.826155e-01     5.289488e-03
    1.5577522687330351e+09     2     6.686786e-01     5.073548e-03
    1.5577524838208699e+09     3     6.422890e-01     4.696714e-03
    1.5577526989721639e+09     4     5.961951e-01     3.871112e-03
    1.5577529139031551e+09     5     5.304139e-01     2.181979e-03
    1.5577531291157911e+09     6     5.000000e-01     5.799814e-06
    1.5577531482987349e+09     7     4.999999e-01     8.587957e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "FD1"
time_iter_q_normgradq = [
    1.557755751472288e9        -1         NaN             NaN
    1.5577561864512529e+09     0     6.898228e-01     5.399008e-03
    1.5577565166688740e+09     1     6.826155e-01     5.289482e-03
    1.5577567876858399e+09     2     6.686786e-01     5.073548e-03
    1.5577570475339010e+09     3     6.422890e-01     4.696713e-03
    1.5577573065499029e+09     4     5.961951e-01     3.871111e-03
    1.5577575670741940e+09     5     5.304139e-01     2.181978e-03
    1.5577578263448451e+09     6     5.000000e-01     5.799812e-06
    1.5577578455805471e+09     7     4.999999e-01     8.598825e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "CSD"
time_iter_q_normgradq = [
    1.557753148494928e9        -1         NaN             NaN
    1.5577535299909971e+09     0     6.898228e-01     5.399008e-03
    1.5577538268994379e+09     1     6.826155e-01     5.289488e-03
    1.5577542159338989e+09     2     6.686786e-01     5.073548e-03
    1.5577545948268349e+09     3     6.422890e-01     4.696714e-03
    1.5577549747877271e+09     4     5.961951e-01     3.871112e-03
    1.5577553535265300e+09     5     5.304139e-01     2.181979e-03
    1.5577557321637819e+09     6     5.000000e-01     5.799814e-06
    1.5577557513279841e+09     7     4.999999e-01     8.587957e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "HYPER"
time_iter_q_normgradq = [
    1.557757845716484e9        -1         NaN             NaN
    1.5577586551774981e+09     0     6.898228e-01     5.399008e-03
    1.5577593583715999e+09     1     6.826155e-01     5.289505e-03
    1.5577603536003621e+09     2     6.686790e-01     5.073625e-03
    1.5577613409561210e+09     3     6.422877e-01     4.696844e-03
    1.5577621722386880e+09     4     5.961760e-01     3.870520e-03
    1.5577629825856421e+09     5     5.304040e-01     2.181599e-03
    1.5577637472446270e+09     6     5.000000e-01     5.804753e-06
    1.5577638844919391e+09     7     4.999999e-01     8.602397e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "FD2"
time_iter_q_normgradq = [
    1.557763884532753e9        -1         NaN             NaN
    1.5577665292074981e+09     0     6.898228e-01     5.399009e-03
    1.5577687054615991e+09     1     6.826159e-01     5.288241e-03
    1.5577697548073909e+09     2     6.687054e-01     5.074894e-03
    1.5577707208770411e+09     3     6.422865e-01     4.698094e-03
    1.5577716878092561e+09     4     5.960494e-01     3.866551e-03
    1.5577726530092020e+09     5     5.303386e-01     2.178785e-03
    1.5577736179376190e+09     6     5.000000e-01     5.858702e-06
    1.5577737469670100e+09     7     4.999999e-01     6.364353e-10
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])


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
lightest_gray   = 0.93 * [1, 1, 1]
lighter_gray   = 0.85 * [1, 1, 1]
darker_gray    = 0.75 * [1, 1, 1]


using Colors
# function to transform RGB (from 0 to 1) color array into hex string
function convert_rgb_to_hex(aRGB)
    r, g, b = aRGB
    return "#" * hex(RGB(r, g, b))
end
colAD2 = convert_rgb_to_hex(lighter_purple)
colAD = convert_rgb_to_hex(lighter_green)
colF1 = convert_rgb_to_hex(darker_red)
colGrid = convert_rgb_to_hex(lightest_gray)


#===========================================
Plot 
===========================================#
fig = figure("pyplot_step_plot", figsize=(8,5))
list_methods = ["F-1", "DUAL", "FD1", "CSD", "HYPER", "FD2" ]
markers      = ["*"  , "D"   , "s"  , "o"  , "D"    , "s"   ]
markersizes  = [9    , 4     , 4    , 4    , 4      , 4     ]
colors       = [colF1, colAD , colAD, colAD, colAD2 , colAD2]
p = Dict()
for (m, mk, ms, c) in zip(list_methods, markers, markersizes, colors)
    #pm = PyPlot.step(x[m], y[m], where="post", marker=mk, markersize=ms, color=c, markeredgecolor="k")
    pm = PyPlot.step(x[m], y[m], where="post", color=c)
    pmend = PyPlot.step(x[m][end], y[m][end], where="post", marker=mk, markersize=ms, color=c, markeredgecolor="k")
    push!(p, m => pm) # needed for legend
    text(x[m][end], y[m][end]/3, m, rotation=90, horizontalalignment="center")
end
#loglog()
semilogy()
xlabel("Computation time (seconds)")
ylabel("Gradient norm")
yticks(10.0 .^ collect(-15:5))
ylim((2e-12, 1e-2))
#xlim((1e2, 1.1e4))
xlim((0, 1.1e4))
#legend(([p[m][1] for m in list_methods]...,), (list_methods...,))
grid("on", which="major", axis="y", color=colGrid)
grid("on", which="minor", axis="x", color=colGrid)
grid("on", which="major", axis="x", color=colGrid)

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
eps_file = joinpath(path_to_package_root, "fig", "Optim_callback.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "Optim_callback.png")
savefig(png_file)


#=

p = df |>
@vlplot(
    width=500,
    height=400,
    mark={
        :line,
        interpolate="step-before",
    },
    encoding={
        x={:time, title="Elapsed computation time (seconds)", scale={typ=:log,domain=[100,30000]}},
        y={:normgradq, title="|∇f|", scale={typ=:log, domain=[1e-11,1e-2]}},
        shape={:method, typ="nominal", scale={range=mymarkers, domain=list_methods}, legend={orient="top-right"}},
        color={:method, typ="nominal", scale={range=mycolors, domain=list_methods}}
    },
    resolve={
        scale={
            color="independent",
            shape="independent"
        },
    }
)

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
pdf_file = joinpath(path_to_package_root, "fig", "Optim_callback_katana_vegalite.pdf")
save(pdf_file, p)


svg_file = joinpath(path_to_package_root, "fig", "Optim_callback_katana_vegalite.svg")
save(svg_file, p)
#=

copied and cleaned from Katana /cluster_output

┌────────────────────────
│ Optimizing using method AF1, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557750778772536e9        -1         NaN             NaN
│    1.5577509205414901e+09     0     6.898228e-01     5.399008e-03
│    1.5577509847967460e+09     1     6.826155e-01     5.289488e-03
│    1.5577510282755060e+09     2     6.686786e-01     5.073548e-03
│    1.5577510610712330e+09     3     6.422890e-01     4.696714e-03
│    1.5577510939584720e+09     4     5.961951e-01     3.871112e-03
│    1.5577511265238740e+09     5     5.304139e-01     2.181979e-03
│    1.5577511590961030e+09     6     5.000000e-01     5.799814e-06
│    1.5577511803307281e+09     7     4.999999e-01     8.587957e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method F1, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557751180373945e9        -1         NaN             NaN
│    1.5577513210150471e+09     0     6.898228e-01     5.399008e-03
│    1.5577513850427470e+09     1     6.826155e-01     5.289488e-03
│    1.5577514285317619e+09     2     6.686786e-01     5.073548e-03
│    1.5577514612624371e+09     3     6.422890e-01     4.696714e-03
│    1.5577514939592800e+09     4     5.961951e-01     3.871112e-03
│    1.5577515265045810e+09     5     5.304139e-01     2.181979e-03
│    1.5577515590177331e+09     6     5.000000e-01     5.799814e-06
│    1.5577515800671041e+09     7     4.999999e-01     8.587957e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method DUAL, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557751580222345e9        -1         NaN             NaN
│    1.5577518531085820e+09     0     6.898228e-01     5.399008e-03
│    1.5577520433873570e+09     1     6.826155e-01     5.289488e-03
│    1.5577522687330351e+09     2     6.686786e-01     5.073548e-03
│    1.5577524838208699e+09     3     6.422890e-01     4.696714e-03
│    1.5577526989721639e+09     4     5.961951e-01     3.871112e-03
│    1.5577529139031551e+09     5     5.304139e-01     2.181979e-03
│    1.5577531291157911e+09     6     5.000000e-01     5.799814e-06
│    1.5577531482987349e+09     7     4.999999e-01     8.587957e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method CSD, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557753148494928e9        -1         NaN             NaN
│    1.5577535299909971e+09     0     6.898228e-01     5.399008e-03
│    1.5577538268994379e+09     1     6.826155e-01     5.289488e-03
│    1.5577542159338989e+09     2     6.686786e-01     5.073548e-03
│    1.5577545948268349e+09     3     6.422890e-01     4.696714e-03
│    1.5577549747877271e+09     4     5.961951e-01     3.871112e-03
│    1.5577553535265300e+09     5     5.304139e-01     2.181979e-03
│    1.5577557321637819e+09     6     5.000000e-01     5.799814e-06
│    1.5577557513279841e+09     7     4.999999e-01     8.587957e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method FD1, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557755751472288e9        -1         NaN             NaN
│    1.5577561864512529e+09     0     6.898228e-01     5.399008e-03
│    1.5577565166688740e+09     1     6.826155e-01     5.289482e-03
│    1.5577567876858399e+09     2     6.686786e-01     5.073548e-03
│    1.5577570475339010e+09     3     6.422890e-01     4.696713e-03
│    1.5577573065499029e+09     4     5.961951e-01     3.871111e-03
│    1.5577575670741940e+09     5     5.304139e-01     2.181978e-03
│    1.5577578263448451e+09     6     5.000000e-01     5.799812e-06
│    1.5577578455805471e+09     7     4.999999e-01     8.598825e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method HYPER, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557757845716484e9        -1         NaN             NaN
│    1.5577586551774981e+09     0     6.898228e-01     5.399008e-03
│    1.5577593583715999e+09     1     6.826155e-01     5.289505e-03
│    1.5577603536003621e+09     2     6.686790e-01     5.073625e-03
│    1.5577613409561210e+09     3     6.422877e-01     4.696844e-03
│    1.5577621722386880e+09     4     5.961760e-01     3.870520e-03
│    1.5577629825856421e+09     5     5.304040e-01     2.181599e-03
│    1.5577637472446270e+09     6     5.000000e-01     5.804753e-06
│    1.5577638844919391e+09     7     4.999999e-01     8.602397e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method FD2, for Precompiled run
│
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.557763884532753e9        -1         NaN             NaN
│    1.5577665292074981e+09     0     6.898228e-01     5.399009e-03
│    1.5577687054615991e+09     1     6.826159e-01     5.288241e-03
│    1.5577697548073909e+09     2     6.687054e-01     5.074894e-03
│    1.5577707208770411e+09     3     6.422865e-01     4.698094e-03
│    1.5577716878092561e+09     4     5.960494e-01     3.866551e-03
│    1.5577726530092020e+09     5     5.303386e-01     2.178785e-03
│    1.5577736179376190e+09     6     5.000000e-01     5.858702e-06
│    1.5577737469670100e+09     7     4.999999e-01     6.364353e-10
└────────────────────────=#

=#
