using VegaLite, DataFrames #
using LaTeXStrings




df = DataFrame(
    method = Array{String}(undef, 0),
    time = Array{Float64}(undef, 0),
    iteration = Array{Int64}(undef, 0),
    qval = Array{Float64}(undef, 0),
    normgradq = Array{Float64}(undef, 0)
)
function mypush!(df, method_name, time_iter_q_normgradq)
    start_time = time_iter_q_normgradq[1, 1]
    for i in 1:size(time_iter_q_normgradq, 1)
        time_iter_q_normgradq[i, 1] -= start_time
        push!(df, [method_name; time_iter_q_normgradq[i, :]])
    end
    return df
end

println("The data being plotted was copy-pasted from Katana output!")
#    time                   iteration      q(λ)           |Dq(λ)|
method_name = "FLASH"
time_iter_q_normgradq = [
    1.552435125178881e9     -1     9.973847e-02     9.846679e-02
    1552435225.5530190468     0     9.973847e-02     9.846679e-02
    1552435324.6035950184     1     1.679872e-02     6.089963e-02
    1552435450.8082480431     2     6.970948e-03     2.493615e-02
    1552435506.2061679363     3     5.678646e-03     2.274644e-03
    1552435551.9947800636     4     5.529870e-03     9.008393e-04
    1552435605.9743819237     5     5.521264e-03     4.901419e-05
    1552435646.1014139652     6     5.521236e-03     5.834945e-07
    1552435673.1531488895     7     5.521236e-03     3.852644e-10
]
mypush!(df, method_name, time_iter_q_normgradq)


method_name = "FiniteDiff"
time_iter_q_normgradq = [
    1.552435687005845e9     -1     9.973847e-02     9.846679e-02
    1552436051.2689580917     0     9.973847e-02     9.846679e-02
    1552436411.1196889877     1     1.679853e-02     6.089921e-02
    1552436817.3611578941     2     6.970752e-03     2.493435e-02
    1552437141.5063579082     3     5.678668e-03     2.275158e-03
    1552437450.2771229744     4     5.529873e-03     9.009779e-04
    1552437772.1610040665     5     5.521264e-03     4.902208e-05
    1552438083.5772459507     6     5.521236e-03     5.855503e-07
    1552438111.2665600777     7     5.521236e-03     3.754357e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "Dual"
time_iter_q_normgradq = [
    1.552438119774152e9     -1     9.973847e-02     9.846679e-02
    1552438408.6534569263     0     9.973847e-02     9.846679e-02
    1552438650.4196999073     1     1.679872e-02     6.089963e-02
    1552438974.5296669006     2     6.970948e-03     2.493615e-02
    1552439175.9353981018     3     5.678646e-03     2.274644e-03
    1552439365.1199560165     4     5.529870e-03     9.008393e-04
    1552439612.5353479385     5     5.521264e-03     4.901419e-05
    1552439800.2384989262     6     5.521236e-03     5.834920e-07
    1552439827.3368079662     7     5.521236e-03     3.873180e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "Complex"
time_iter_q_normgradq = [
    1.552439835793842e9     -1     9.973847e-02     9.846679e-02
    1552440299.3248779774     0     9.973847e-02     9.846679e-02
    1552440674.7697989941     1     1.679872e-02     6.089963e-02
    1552441186.7565340996     2     6.970948e-03     2.493615e-02
    1552441514.2213349342     3     5.678646e-03     2.274644e-03
    1552441837.0221610069     4     5.529870e-03     9.008393e-04
    1552442275.9999530315     5     5.521264e-03     4.901418e-05
    1552442603.1407780647     6     5.521236e-03     5.834949e-07
    1552442630.4078490734     7     5.521236e-03     3.976156e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "HyperDual"
time_iter_q_normgradq = [
    1.552442639194264e9     -1     9.973847e-02     9.846679e-02
    1552443115.0007998943     0     9.973847e-02     9.846679e-02
    1552443541.8403270245     1     1.679872e-02     6.089963e-02
    1552444049.6967279911     2     6.970948e-03     2.493615e-02
    1552444433.5782051086     3     5.678646e-03     2.274644e-03
    1552444807.2067389488     4     5.529870e-03     9.008393e-04
    1552445240.5581049919     5     5.521264e-03     4.901419e-05
    1552445612.3695719242     6     5.521236e-03     5.834910e-07
    1552445732.8631589413     7     5.521236e-03     3.690092e-10
]
mypush!(df, method_name, time_iter_q_normgradq)


method_name = "sym-HyperDual"
time_iter_q_normgradq = [
    1.552445741551325e9     -1     9.973847e-02     9.846679e-02
    1552446100.4548881054     0     9.973847e-02     9.846679e-02
    1552446419.1283090115     1     1.679872e-02     6.089963e-02
    1552446819.625041008     2     6.970948e-03     2.493615e-02
    1552447097.2709009647     3     5.678646e-03     2.274644e-03
    1552447362.8477189541     4     5.529870e-03     9.008393e-04
    1552447686.9819951057     5     5.521264e-03     4.901419e-05
    1552447951.0648400784     6     5.521236e-03     5.834910e-07
    1552448068.9834229946     7     5.521236e-03     3.690092e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "buf-HyperDual"
time_iter_q_normgradq = [
    1.552448077731379e9     -1     9.973847e-02     9.846679e-02
    1552448292.5661931038     0     9.973847e-02     9.846679e-02
    1552448515.7086060047     1     1.679872e-02     6.089963e-02
    1552448767.844135046     2     6.970948e-03     2.493615e-02
    1552448946.2576050758     3     5.678646e-03     2.274644e-03
    1552449115.9403719902     4     5.529870e-03     9.008393e-04
    1552449294.626003027     5     5.521264e-03     4.901418e-05
    1552449468.7211060524     6     5.521236e-03     5.834946e-07
    1552449543.2695729733     7     5.521236e-03     3.852952e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

list_methods = ["FLASH", "Dual", "buf-HyperDual", "sym-HyperDual", "HyperDual", "FiniteDiff", "Complex"]

Array_colors_RGB = [
     0   0   0
     0 114 178
    66 160 213
    86 180 233
   106 200 253
     0 158 115
   230 159   0
]

using Colors
# function to transform RGB (from 0 to 255) color array into list of hex
function Array_color_RGB_to_hex_list(M)
    out = Array{String,1}(undef, 0)
    for icol in 1:size(M, 1)
        r, g, b = M[icol, :] / 255
        push!(out, "#" * hex(RGB(r, g, b)))
    end 
    return out
end

mycolors = Array_color_RGB_to_hex_list(Array_colors_RGB)

p = df |>
@vlplot(
    width=300,
    mark={
        :line,
        interpolate="step-before"
    },
    encoding={
        x={:time, title="Elapsed computation time (seconds)"},
        y={:normgradq, title="norm of gradient", scale={typ=:log}},
        shape={:method, typ="nominal", scale={domain=list_methods}},
        color={:method, typ="nominal", scale={range=mycolors, domain=list_methods}}
    },
    resolve={
        scale={
            color="independent",
            shape="independent"
        },
    }
)

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
pdf_file = joinpath(path_to_package_root, "fig", "Optim_callback_katana_vegalite.pdf")
save(pdf_file, p)

#=

┌────────────────────────
│ Optimizing using method FLASH, and printing time and q :
│
│    1.552285149618742e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552285262.5103099346     0     9.973847e-02     9.846679e-02
│    1552285371.4692571163     1     1.679872e-02     6.089963e-02
│    1552285511.096765995     2     6.970948e-03     2.493615e-02
│    1552285575.2337169647     3     5.678646e-03     2.274644e-03
│    1552285627.0646879673     4     5.529870e-03     9.008393e-04
│    1552285683.9633190632     5     5.521264e-03     4.901419e-05
│    1552285726.5941050053     6     5.521236e-03     5.834928e-07
│    1552285755.4307689667     7     5.521236e-03     3.744009e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method FiniteDiff, and printing time and q :
│
│    1.552285769695254e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552286163.6655740738     0     9.973847e-02     9.846679e-02
│    1552286547.5203258991     1     1.679853e-02     6.089921e-02
│    1552286976.8766720295     2     6.970752e-03     2.493435e-02
│    1552287323.8391370773     3     5.678668e-03     2.275158e-03
│    1552287654.4600300789     4     5.529873e-03     9.009779e-04
│    1552287999.3380770683     5     5.521264e-03     4.902208e-05
│    1552288326.3779881001     6     5.521236e-03     5.855503e-07
│    1552288357.3542480469     7     5.521236e-03     3.754357e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method Dual, and printing time and q :
│
│    1.552288367126222e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552288687.195045948     0     9.973847e-02     9.846679e-02
│    1552288954.4081599712     1     1.679872e-02     6.089963e-02
│    1552289317.5387248993     2     6.970948e-03     2.493615e-02
│    1552289539.5822749138     3     5.678646e-03     2.274644e-03
│    1552289757.4462211132     4     5.529870e-03     9.008393e-04
│    1552290089.9508259296     5     5.521264e-03     4.901418e-05
│    1552290305.0491220951     6     5.521236e-03     5.834965e-07
│    1552290333.1003849506     7     5.521236e-03     3.862567e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method Complex, and printing time and q :
│
│    1.552290341901904e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552290847.9102509022     0     9.973847e-02     9.846679e-02
│    1552291250.8377010822     1     1.679872e-02     6.089963e-02
│    1552291806.8812630177     2     6.970948e-03     2.493615e-02
│    1552292157.0304598808     3     5.678646e-03     2.274644e-03
│    1552292503.4707729816     4     5.529870e-03     9.008393e-04
│    1552292975.3528540134     5     5.521264e-03     4.901420e-05
│    1552293329.6839270592     6     5.521236e-03     5.834909e-07
│    1552293357.8559169769     7     5.521236e-03     3.617950e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method HyperDual, and printing time and q :
│
│    1.552293366675875e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552294654.9384651184     0     9.973847e-02     9.846679e-02
│    1552295606.6262910366     1     1.679872e-02     6.089963e-02
│    1552296939.1215419769     2     6.970948e-03     2.493615e-02
│    1552297865.3929319382     3     5.678646e-03     2.274644e-03
│    1552298757.0831110477     4     5.529870e-03     9.008393e-04
│    1552300001.1723589897     5     5.521264e-03     4.901419e-05
│    1552300914.7481648922     6     5.521236e-03     5.834912e-07
│    1552301041.0953009129     7     5.521236e-03     3.724027e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method sHyperDual, and printing time and q :
│
│    1.552301050019803e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552301918.547385931     0     9.973847e-02     9.846679e-02
│    1552302578.5258059502     1     1.679872e-02     6.089963e-02
│    1552303495.229377985     2     6.970948e-03     2.493615e-02
│    1552304112.7929279804     3     5.678646e-03     2.274644e-03
│    1552304697.6630148888     4     5.529870e-03     9.008393e-04
│    1552305537.2959148884     5     5.521264e-03     4.901419e-05
│    1552306143.3304190636     6     5.521236e-03     5.834912e-07
│    1552306276.5732889175     7     5.521236e-03     3.724027e-10
└────────────────────────

┌────────────────────────
│ Optimizing using method bHyperDual, and printing time and q :
│
│    1.552306286265865e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1552306997.5528509617     0     9.973847e-02     9.846679e-02
│    1552307552.6609730721     1     1.679872e-02     6.089963e-02
│    1552308528.7847630978     2     6.970948e-03     2.493615e-02
│    1552309026.8967139721     3     5.678646e-03     2.274644e-03
│    1552309499.8328630924     4     5.529870e-03     9.008393e-04
│    1552310336.5372889042     5     5.521264e-03     4.901418e-05
│    1552310812.5246040821     6     5.521236e-03     5.834950e-07
│    1552310988.0982470512     7     5.521236e-03     3.830555e-10
└────────────────────────

=#
