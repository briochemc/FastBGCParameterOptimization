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
    1.558026139908681e9        -1         NaN            NaN
    1.5580261992445560e+09     0     8.450480e-02     3.445963e-01
    1.5580262540682299e+09     1     1.093143e-02     2.234148e-02
    1.5580262909548481e+09     2     1.093143e-02     2.234148e-02
    1.5580263362981789e+09     3     7.361795e-03     1.693728e-02
    1.5580264004758489e+09     4     5.848101e-03     4.571342e-03
    1.5580264561298239e+09     5     4.038817e-03     2.709566e-03
    1.5580264949480760e+09     6     2.577054e-03     2.776576e-03
    1.5580265297423429e+09     7     2.533827e-03     4.516446e-05
    1.5580265641108410e+09     8     2.533813e-03     2.173716e-07
    1.5580265863829119e+09     9     2.533813e-03     4.962021e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "DUAL"
time_iter_q_normgradq = [
    1.558026586389074e9        -1         NaN            NaN
    1.5580267725074501e+09     0     8.450480e-02     3.445963e-01
    1.5580269588302510e+09     1     1.093143e-02     2.234148e-02
    1.5580269940333941e+09     2     1.093143e-02     2.234148e-02
    1.5580272232405550e+09     3     7.361795e-03     1.693728e-02
    1.5580274755248849e+09     4     5.848101e-03     4.571342e-03
    1.5580276608299589e+09     5     4.038817e-03     2.709566e-03
    1.5580278276107180e+09     6     2.577054e-03     2.776576e-03
    1.5580280493897791e+09     7     2.533827e-03     4.516447e-05
    1.5580282713219290e+09     8     2.533813e-03     2.172284e-07
    1.5580282916157739e+09     9     2.533813e-03     8.054451e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "FD1"
time_iter_q_normgradq = [
    1.558031118459424e9        -1         NaN            NaN
    1.5580314191674490e+09     0     8.450480e-02     3.445963e-01
    1.5580317158880091e+09     1     1.093143e-02     2.234148e-02
    1.5580317504852829e+09     2     1.093143e-02     2.234148e-02
    1.5580320342729540e+09     3     7.361795e-03     1.693727e-02
    1.5580323368482330e+09     4     5.848101e-03     4.571342e-03
    1.5580326328577459e+09     5     4.038820e-03     2.709413e-03
    1.5580329112136080e+09     6     2.577055e-03     2.776572e-03
    1.5580331849400189e+09     7     2.533827e-03     4.519602e-05
    1.5580334572237720e+09     8     2.533813e-03     2.208897e-07
    1.5580334773777430e+09     9     2.533813e-03     1.610760e-10
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "CSD"
time_iter_q_normgradq = [
    1.558028291743703e9        -1         NaN            NaN
    1.5580285938239219e+09     0     8.450480e-02     3.445963e-01
    1.5580288897053709e+09     1     1.093143e-02     2.234148e-02
    1.5580289237822220e+09     2     1.093143e-02     2.234148e-02
    1.5580293256107790e+09     3     7.361795e-03     1.693728e-02
    1.5580297442267101e+09     4     5.848101e-03     4.571342e-03
    1.5580300433359160e+09     5     4.038817e-03     2.709566e-03
    1.5580303206879370e+09     6     2.577054e-03     2.776576e-03
    1.5580307089946721e+09     7     2.533827e-03     4.516445e-05
    1.5580310981021700e+09     8     2.533813e-03     2.173728e-07
    1.5580311183136220e+09     9     2.533813e-03     9.758266e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "HYPER"
time_iter_q_normgradq = [
    1.558033477834328e9        -1         NaN            NaN
    1.5580341928867800e+09     0     8.450480e-02     3.445963e-01
    1.5580349044803801e+09     1     1.093143e-02     2.234148e-02
    1.5580350032748921e+09     2     1.093143e-02     2.234148e-02
    1.5580359923737810e+09     3     7.361793e-03     1.693714e-02
    1.5580370365108819e+09     4     5.848101e-03     4.571341e-03
    1.5580377409823639e+09     5     4.038817e-03     2.709567e-03
    1.5580384329587841e+09     6     2.577054e-03     2.776576e-03
    1.5580394468175941e+09     7     2.533827e-03     4.516757e-05
    1.5580404548681681e+09     8     2.533813e-03     2.173940e-07
    1.5580405933614130e+09     9     2.533813e-03     6.190045e-11
]
push!(x, m => time_iter_q_normgradq[2:end,1] .- time_iter_q_normgradq[1,1])
push!(y, m => time_iter_q_normgradq[2:end,4])

m = "FD2"
time_iter_q_normgradq = [
    1.558040593444234e9        -1         NaN            NaN
    1.5580417865699401e+09     0     8.450480e-02     3.445964e-01
    1.5580429702467480e+09     1     1.093147e-02     2.234103e-02
    1.5580431215995660e+09     2     1.093147e-02     2.234103e-02
    1.5580443344388211e+09     3     7.361233e-03     1.692947e-02
    1.5580455642303290e+09     4     5.847841e-03     4.570866e-03
    1.5580467310326970e+09     5     4.038679e-03     2.708307e-03
    1.5580477861093731e+09     6     2.576974e-03     2.776960e-03
    1.5580488365545111e+09     7     2.533862e-03     4.471783e-05
    1.5580498831764319e+09     8     2.533814e-03     6.991567e-07
    1.5580509299039240e+09     9     2.533813e-03     3.269161e-08
    1.5580510647831969e+09    10     2.533813e-03     3.321707e-09
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
    text(x[m][end], y[m][end]/4, m, rotation=90, horizontalalignment="center")
end
#loglog()
semilogy()
xlabel("Computation time (seconds)")
ylabel("Gradient norm")
yticks(10.0 .^ collect(-15:5))
ylim((1e-12, 1e0))
#xlim((1e2, 1.1e4))
xlim((0, 1.1e4))
#legend(([p[m][1] for m in list_methods]...,), (list_methods...,))
grid("on", which="major", axis="y", color=colGrid)
grid("on", which="minor", axis="x", color=colGrid)
grid("on", which="major", axis="x", color=colGrid)
plot([xlim()...], 1e-8 * [1, 1], color=darker_gray)

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
eps_file = joinpath(path_to_package_root, "fig", "Optim_callback.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "Optim_callback.png")
savefig(png_file)


#=

copied and cleaned from Katana /cluster_output

┌────────────────────────
│ Optimizing using method AF1, for Precompiled run
│
│    1.558025645358724e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580257123228760e+09     0     8.450480e-02     3.445963e-01
│    1.5580257727929730e+09     1     1.093143e-02     2.234148e-02
│    1.5580258121005919e+09     2     1.093143e-02     2.234148e-02
│    1.5580258640124340e+09     3     7.361795e-03     1.693728e-02
│    1.5580259320926700e+09     4     5.848101e-03     4.571342e-03
│    1.5580259944940870e+09     5     4.038817e-03     2.709566e-03
│    1.5580260388522351e+09     6     2.577054e-03     2.776576e-03
│    1.5580260783906600e+09     7     2.533827e-03     4.516430e-05
│    1.5580261170834539e+09     8     2.533813e-03     2.174209e-07
│    1.5580261398953979e+09     9     2.533813e-03     3.265780e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method F1, for Precompiled run
│
│    1.558026139908681e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580261992445560e+09     0     8.450480e-02     3.445963e-01
│    1.5580262540682299e+09     1     1.093143e-02     2.234148e-02
│    1.5580262909548481e+09     2     1.093143e-02     2.234148e-02
│    1.5580263362981789e+09     3     7.361795e-03     1.693728e-02
│    1.5580264004758489e+09     4     5.848101e-03     4.571342e-03
│    1.5580264561298239e+09     5     4.038817e-03     2.709566e-03
│    1.5580264949480760e+09     6     2.577054e-03     2.776576e-03
│    1.5580265297423429e+09     7     2.533827e-03     4.516446e-05
│    1.5580265641108410e+09     8     2.533813e-03     2.173716e-07
│    1.5580265863829119e+09     9     2.533813e-03     4.962021e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method DUAL, for Precompiled run
│
│    1.558026586389074e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580267725074501e+09     0     8.450480e-02     3.445963e-01
│    1.5580269588302510e+09     1     1.093143e-02     2.234148e-02
│    1.5580269940333941e+09     2     1.093143e-02     2.234148e-02
│    1.5580272232405550e+09     3     7.361795e-03     1.693728e-02
│    1.5580274755248849e+09     4     5.848101e-03     4.571342e-03
│    1.5580276608299589e+09     5     4.038817e-03     2.709566e-03
│    1.5580278276107180e+09     6     2.577054e-03     2.776576e-03
│    1.5580280493897791e+09     7     2.533827e-03     4.516447e-05
│    1.5580282713219290e+09     8     2.533813e-03     2.172284e-07
│    1.5580282916157739e+09     9     2.533813e-03     8.054451e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method CSD, for Precompiled run
│
│    1.558028291743703e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580285938239219e+09     0     8.450480e-02     3.445963e-01
│    1.5580288897053709e+09     1     1.093143e-02     2.234148e-02
│    1.5580289237822220e+09     2     1.093143e-02     2.234148e-02
│    1.5580293256107790e+09     3     7.361795e-03     1.693728e-02
│    1.5580297442267101e+09     4     5.848101e-03     4.571342e-03
│    1.5580300433359160e+09     5     4.038817e-03     2.709566e-03
│    1.5580303206879370e+09     6     2.577054e-03     2.776576e-03
│    1.5580307089946721e+09     7     2.533827e-03     4.516445e-05
│    1.5580310981021700e+09     8     2.533813e-03     2.173728e-07
│    1.5580311183136220e+09     9     2.533813e-03     9.758266e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method FD1, for Precompiled run
│
│    1.558031118459424e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580314191674490e+09     0     8.450480e-02     3.445963e-01
│    1.5580317158880091e+09     1     1.093143e-02     2.234148e-02
│    1.5580317504852829e+09     2     1.093143e-02     2.234148e-02
│    1.5580320342729540e+09     3     7.361795e-03     1.693727e-02
│    1.5580323368482330e+09     4     5.848101e-03     4.571342e-03
│    1.5580326328577459e+09     5     4.038820e-03     2.709413e-03
│    1.5580329112136080e+09     6     2.577055e-03     2.776572e-03
│    1.5580331849400189e+09     7     2.533827e-03     4.519602e-05
│    1.5580334572237720e+09     8     2.533813e-03     2.208897e-07
│    1.5580334773777430e+09     9     2.533813e-03     1.610760e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method HYPER, for Precompiled run
│
│    1.558033477834328e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580341928867800e+09     0     8.450480e-02     3.445963e-01
│    1.5580349044803801e+09     1     1.093143e-02     2.234148e-02
│    1.5580350032748921e+09     2     1.093143e-02     2.234148e-02
│    1.5580359923737810e+09     3     7.361793e-03     1.693714e-02
│    1.5580370365108819e+09     4     5.848101e-03     4.571341e-03
│    1.5580377409823639e+09     5     4.038817e-03     2.709567e-03
│    1.5580384329587841e+09     6     2.577054e-03     2.776576e-03
│    1.5580394468175941e+09     7     2.533827e-03     4.516757e-05
│    1.5580404548681681e+09     8     2.533813e-03     2.173940e-07
│    1.5580405933614130e+09     9     2.533813e-03     6.190045e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method FD2, for Precompiled run
│
│    1.558040593444234e9
│    time                   iteration      f̂(λ)           |∇f̂(λ)|
│    1.5580417865699401e+09     0     8.450480e-02     3.445964e-01
│    1.5580429702467480e+09     1     1.093147e-02     2.234103e-02
│    1.5580431215995660e+09     2     1.093147e-02     2.234103e-02
│    1.5580443344388211e+09     3     7.361233e-03     1.692947e-02
│    1.5580455642303290e+09     4     5.847841e-03     4.570866e-03
│    1.5580467310326970e+09     5     4.038679e-03     2.708307e-03
│    1.5580477861093731e+09     6     2.576974e-03     2.776960e-03
│    1.5580488365545111e+09     7     2.533862e-03     4.471783e-05
│    1.5580498831764319e+09     8     2.533814e-03     6.991567e-07
│    1.5580509299039240e+09     9     2.533813e-03     3.269161e-08
│    1.5580510647831969e+09    10     2.533813e-03     3.321707e-09
└────────────────────────

=#
