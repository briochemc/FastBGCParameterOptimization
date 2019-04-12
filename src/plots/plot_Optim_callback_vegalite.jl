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
    for i in 2:size(time_iter_q_normgradq, 1)
        time_iter_q_normgradq[i, 1] -= start_time
        time_iter_q_normgradq[:, 1] .= max.(time_iter_q_normgradq[:, 1], 100.0)
        push!(df, [method_name; time_iter_q_normgradq[i, :]])
    end
    return df
end


println("The data being plotted was copy-pasted from Katana output!")
#    time                   iteration      q(λ)           |Dq(λ)|
method_name = "F-1"
time_iter_q_normgradq = [
    1.554931957502695e9        -1            NaN             NaN
    1.5549320565900841e+09     0     9.973847e-02     9.846679e-02
    1.5549321534185779e+09     1     1.679872e-02     6.089963e-02
    1.5549322792901530e+09     2     6.970948e-03     2.493615e-02
    1.5549323296461041e+09     3     5.678646e-03     2.274644e-03
    1.5549323700646110e+09     4     5.529870e-03     9.008393e-04
    1.5549324199735501e+09     5     5.521264e-03     4.901419e-05
    1.5549324544249890e+09     6     5.521236e-03     5.834978e-07
    1.5549324764175520e+09     7     5.521236e-03     3.886968e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "DUAL"
time_iter_q_normgradq = [
    1.554932476422953e9        -1            NaN             NaN
    1.5549327393173470e+09     0     9.973847e-02     9.846679e-02
    1.5549329462172129e+09     1     1.679872e-02     6.089963e-02
    1.5549332329707501e+09     2     6.970948e-03     2.493615e-02
    1.5549333906547461e+09     3     5.678646e-03     2.274644e-03
    1.5549335368875189e+09     4     5.529870e-03     9.008393e-04
    1.5549337428765121e+09     5     5.521264e-03     4.901419e-05
    1.5549338860186059e+09     6     5.521236e-03     5.834973e-07
    1.5549339057977159e+09     7     5.521236e-03     3.841849e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "FD1"
time_iter_q_normgradq = [
    1.554936086725043e9        -1            NaN             NaN
    1.5549363799715650e+09     0     9.973847e-02     9.846679e-02
    1.5549366622112429e+09     1     1.679872e-02     6.089962e-02
    1.5549369976018641e+09     2     6.970946e-03     2.493613e-02
    1.5549372475325570e+09     3     5.678664e-03     2.250981e-03
    1.5549374905341091e+09     4     5.529873e-03     9.013453e-04
    1.5549377440755579e+09     5     5.521264e-03     4.905180e-05
    1.5549379851920121e+09     6     5.521236e-03     5.834383e-07
    1.5549380056898539e+09     7     5.521236e-03     3.792707e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "CSD"
time_iter_q_normgradq = [
    1.554933905967612e9        -1            NaN             NaN
    1.5549342912626989e+09     0     9.973847e-02     9.846679e-02
    1.5549345805700810e+09     1     1.679872e-02     6.089963e-02
    1.5549350051208200e+09     2     6.970948e-03     2.493615e-02
    1.5549352478671119e+09     3     5.678646e-03     2.274644e-03
    1.5549354833078370e+09     4     5.529870e-03     9.008393e-04
    1.5549358301850841e+09     5     5.521264e-03     4.901419e-05
    1.5549360668987601e+09     6     5.521236e-03     5.834950e-07
    1.5549360865820920e+09     7     5.521236e-03     3.858098e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "HYPER"
time_iter_q_normgradq = [
    1.554938005827915e9        -1            NaN             NaN
    1.5549388540773749e+09     0     9.973847e-02     9.846679e-02
    1.5549394392444971e+09     1     1.679872e-02     6.089963e-02
    1.5549403306286790e+09     2     6.970948e-03     2.493615e-02
    1.5549408725899410e+09     3     5.678646e-03     2.274644e-03
    1.5549414042682891e+09     4     5.529870e-03     9.008393e-04
    1.5549422032390299e+09     5     5.521264e-03     4.901418e-05
    1.5549427448844261e+09     6     5.521236e-03     5.834916e-07
    1.5549428184643610e+09     7     5.521236e-03     3.860844e-10
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "FD2"
time_iter_q_normgradq = [
    1.554942818498172e9        -1            NaN             NaN
    1.5549438858878219e+09     0     9.973847e-02     9.846679e-02
    1.5549447314547229e+09     1     1.675709e-02     6.081409e-02
    1.5549463066818390e+09     2     6.927449e-03     2.450954e-02
    1.5549473549851000e+09     3     5.660650e-03     2.003297e-03
    1.5549486106761160e+09     4     5.527942e-03     9.347555e-04
    1.5549499118954070e+09     5     5.521245e-03     3.489643e-05
    1.5549500563287971e+09     6     5.521245e-03     3.489643e-05
    1.5549502035800760e+09     7     5.521245e-03     3.489643e-05
    1.5549503446897969e+09     8     5.521245e-03     3.489643e-05
    1.5549504660686951e+09     9     5.521245e-03     3.489643e-05
    1.5549518515212281e+09    10     5.521242e-03     2.376590e-06
    1.5549532392969739e+09    11     5.521238e-03     1.122145e-06
    1.5549546077921720e+09    12     5.521236e-03     4.092866e-07
    1.5549547270177040e+09    13     5.521236e-03     4.092866e-07
    1.5549560854869180e+09    14     5.521236e-03     1.068310e-07
    1.5549562052175820e+09    15     5.521236e-03     1.068310e-07
    1.5549575683601999e+09    16     5.521236e-03     4.618704e-08
    1.5549576893180699e+09    17     5.521236e-03     4.618704e-08
    1.5549590547731459e+09    18     5.521236e-03     3.069491e-08
    1.5549604088873370e+09    19     5.521236e-03     2.199276e-08
    1.5549617769524150e+09    20     5.521236e-03     1.017175e-08
    1.5549618965016000e+09    21     5.521236e-03     1.017175e-08
    1.5549620164352100e+09    22     5.521236e-03     6.222843e-09
]
mypush!(df, method_name, time_iter_q_normgradq)


# Julia_colors_RGB 
darker_blue    = [0.251, 0.388, 0.847]
lighter_blue   = [0.4  , 0.51 , 0.878]
darker_purple  = [0.584, 0.345, 0.698]
lighter_purple = [0.667, 0.475, 0.757]
darker_green   = [0.22 , 0.596, 0.149]
lighter_green  = [0.376, 0.678, 0.318]
darker_red     = [0.796, 0.235, 0.2  ]
lighter_red    = [0.835, 0.388, 0.361]

mycolors_RGB = 0.85ones(6,3) # gray base
mycolors_RGB[1, :] .= darker_red
mycolors_RGB[2:4, :] .= 0.65
mymarkers = ["cross", "diamond", "square", "circle", "diamond", "square"]
list_methods = ["F-1", "DUAL", "FD1", "CSD", "HYPER", "FD2"]
filter!(row -> row[:method] ∈ list_methods, df)

using Colors
# function to transform RGB (from 0 to 1) color array into list of hex
function array_of_RGB_to_hex_list(M)
    out = Array{String,1}(undef, 0)
    for icol in 1:size(M, 1)
        r, g, b = M[icol, :]
        push!(out, "#" * hex(RGB(r, g, b)))
    end
    return out
end
#
mycolors = array_of_RGB_to_hex_list(mycolors_RGB)
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
        y={:normgradq, title="|∇f|", scale={typ=:log, domain=[1e-11,1e-1]}},
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

───────────────────────
Optimizing using method OF1, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method F0, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method F1, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method DUAL, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method CSD, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method FD1, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method HYPER, for Precompiled run

───────────────────────

───────────────────────
Optimizing using method FD2, for Precompiled run

─────────────────────


=#
