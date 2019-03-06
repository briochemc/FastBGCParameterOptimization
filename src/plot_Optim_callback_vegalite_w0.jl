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
method_name = "D2q"
time_iter_q_normgradq = [
    1.551757506846656e9    -1     6.521377e-02     5.500858e-02
    1551757571.703897953     0     6.521377e-02     5.500858e-02
    1551757624.4161579609     1     1.367760e-02     2.434447e-02
    1551757671.4375801086     2     7.791760e-03     2.402866e-03
    1551757707.9004130363     3     7.791760e-03     2.402866e-03
    1551757752.3280088902     4     7.238546e-03     1.034544e-03
    1551757784.0264921188     5     7.238546e-03     1.034544e-03
    1551757826.9310441017     6     7.228064e-03     1.213862e-04
    1551757864.5752470493     7     7.223945e-03     7.629348e-04
    1551757900.3933150768     8     7.217872e-03     4.117735e-04
    1551757938.3700749874     9     7.214912e-03     6.281950e-04
    1551757966.9808549881    10     7.214912e-03     6.281950e-04
    1551757994.0706369877    11     7.214912e-03     6.281950e-04
    1551758021.1456990242    12     7.214912e-03     6.281950e-04
    1551758058.6457710266    13     7.213113e-03     1.350311e-05
    1551758094.9337079525    14     7.212779e-03     2.790561e-05
    1551758124.6829230785    15     7.212779e-03     2.790561e-05
    1551758161.0637609959    16     7.212673e-03     7.390252e-06
    1551758188.5657811165    17     7.212673e-03     7.390252e-06
    1551758223.2175629139    18     7.212639e-03     1.560990e-06
    1551758257.7572259903    19     7.212608e-03     7.172192e-06
    1551758290.823802948    20     7.212608e-03     5.154701e-08
    1551758316.8794739246    21     7.212608e-03     5.154701e-08
    1551758342.8025970459    22     7.212608e-03     5.154701e-08
    1551758368.7285950184    23     7.212608e-03     5.154701e-08
    1551758394.6267869473    24     7.212608e-03     5.154701e-08
    1551758420.5754930973    25     7.212608e-03     5.154701e-08
    1551758446.5254330635    26     7.212608e-03     5.154701e-08
    1551758472.454957962    27     7.212608e-03     2.368872e-11
]
mypush!(df, method_name, time_iter_q_normgradq)


method_name = "FDDq"
time_iter_q_normgradq = [
    1.551758485966057e9     -1     6.521377e-02     5.500858e-02
    1551758730.3310580254     0     6.521377e-02     5.500858e-02
    1551758980.3469901085     1     1.388960e-02     2.473551e-02
    1551759231.7757260799     2     7.817013e-03     2.440243e-03
    1551759268.4188189507     3     7.817013e-03     2.440243e-03
    1551759513.465708971     4     7.240150e-03     1.108362e-03
    1551759545.0571510792     5     7.240150e-03     1.108362e-03
    1551759790.3847711086     6     7.228540e-03     1.137964e-04
    1551760027.7207050323     7     7.224389e-03     7.296822e-04
    1551760264.1746320724     8     7.218334e-03     4.886249e-04
    1551760502.5484640598     9     7.214718e-03     4.958740e-04
    1551760740.1306209564    10     7.212923e-03     9.511154e-05
    1551760976.5135378838    11     7.212629e-03     7.209988e-05
    1551761210.6786019802    12     7.212608e-03     1.301734e-07
    1551761236.66918993    13     7.212608e-03     1.301734e-07
    1551761262.6539320946    14     7.212608e-03     1.301734e-07
    1551761288.6356539726    15     7.212608e-03     1.301734e-07
    1551761314.6406729221    16     7.212608e-03     1.301734e-07
    1551761340.6277658939    17     7.212608e-03     1.301734e-07
    1551761366.6241810322    18     7.212608e-03     1.301734e-07
    1551761392.6452329159    19     7.212608e-03     1.301734e-07
    1551761418.6504468918    20     7.212608e-03     1.147283e-08   
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "ADDq"
time_iter_q_normgradq = [
    1.551761427202471e9     -1     6.521377e-02     5.500858e-02
    1551761629.6983449459     0     6.521377e-02     5.500858e-02
    1551761940.2861320972     1     1.367760e-02     2.434447e-02
    1551762274.8026890755     2     7.791760e-03     2.402866e-03
    1551762311.2957570553     3     7.791760e-03     2.402866e-03
    1551762491.3001790047     4     7.238546e-03     1.034544e-03
    1551762531.9407939911     5     7.238546e-03     1.034544e-03
    1551763008.6284229755     6     7.228064e-03     1.213862e-04
    1551763434.4954199791     7     7.223945e-03     7.629348e-04
    1551763759.7934660912     8     7.217872e-03     4.117736e-04
    1551763789.5299289227     9     7.217872e-03     4.117736e-04
    1551764300.1290330887    10     7.214556e-03     4.384884e-04
    1551764511.0417521    11     7.212904e-03     1.137828e-04
    1551764789.6769599915    12     7.212623e-03     5.275763e-05
    1551764815.7157950401    13     7.212623e-03     5.275763e-05
    1551764841.7746610641    14     7.212623e-03     5.275763e-05
    1551764867.7038359642    15     7.212623e-03     5.275763e-05
    1551764893.7761321068    16     7.212623e-03     5.275763e-05
    1551765233.697437048    17     7.212608e-03     8.153595e-07
    1551765504.6984729767    18     7.212608e-03     2.052311e-08
    1551765530.7620489597    19     7.212608e-03     1.092557e-11
]
mypush!(df, method_name, time_iter_q_normgradq)

method_name = "CSDDq"
time_iter_q_normgradq = [
    1.551765539320366e9     -1     6.521377e-02     5.500858e-02
    1551766004.9462089539     0     6.521377e-02     5.500858e-02
    1551766644.1099970341     1     1.367760e-02     2.434447e-02
    1551767213.5281050205     2     7.791760e-03     2.402866e-03
    1551767250.170830965     3     7.791760e-03     2.402866e-03
    1551767950.7097389698     4     7.238546e-03     1.034544e-03
    1551767982.7712020874     5     7.238546e-03     1.034544e-03
    1551768626.3271861076     6     7.228064e-03     1.213862e-04
    1551769053.2136509418     7     7.223945e-03     7.629349e-04
    1551769380.238038063     8     7.217872e-03     4.117735e-04
    1551769410.0014419556     9     7.217872e-03     4.117735e-04
    1551769956.7617139816    10     7.214556e-03     4.384883e-04
    1551769985.4385509491    11     7.214556e-03     4.384883e-04
    1551770404.5565180779    12     7.213443e-03     1.422530e-05
    1551770434.0334420204    13     7.213443e-03     1.422530e-05
    1551771309.8578929901    14     7.212996e-03     2.642094e-05
    1551771816.9819290638    15     7.212656e-03     1.104680e-04
    1551772265.233741045    16     7.212608e-03     1.281638e-08
    1551772291.276966095    17     7.212608e-03     3.170371e-11
]
mypush!(df, method_name, time_iter_q_normgradq)


list_methods = [
    "D2q"
    "FDDq"
    "ADDq"
    "CSDDq"
]

p = df |>
@vlplot(
    mark={
        :line,
        interpolate="step-before"
    },
    x={:time, title="Elapsed computation time (seconds)"},
    y={:normgradq, title="norm of dq(p)", scale={typ=:log}},
    color={:method, scale={scheme="category10"}}
)

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
pdf_file = joinpath(path_to_package_root, "fig", "Optim_callback_w0_katana_vegalite.pdf")
save(pdf_file, p)

#=


┌────────────────────────
│ Optimizing using method D2q, and printing time and q :
│
│    1.551757506846656e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1551757571.703897953     0     6.521377e-02     5.500858e-02
│    1551757624.4161579609     1     1.367760e-02     2.434447e-02
│    1551757671.4375801086     2     7.791760e-03     2.402866e-03
│    1551757707.9004130363     3     7.791760e-03     2.402866e-03
│    1551757752.3280088902     4     7.238546e-03     1.034544e-03
│    1551757784.0264921188     5     7.238546e-03     1.034544e-03
│    1551757826.9310441017     6     7.228064e-03     1.213862e-04
│    1551757864.5752470493     7     7.223945e-03     7.629348e-04
│    1551757900.3933150768     8     7.217872e-03     4.117735e-04
│    1551757938.3700749874     9     7.214912e-03     6.281950e-04
│    1551757966.9808549881    10     7.214912e-03     6.281950e-04
│    1551757994.0706369877    11     7.214912e-03     6.281950e-04
│    1551758021.1456990242    12     7.214912e-03     6.281950e-04
│    1551758058.6457710266    13     7.213113e-03     1.350311e-05
│    1551758094.9337079525    14     7.212779e-03     2.790561e-05
│    1551758124.6829230785    15     7.212779e-03     2.790561e-05
│    1551758161.0637609959    16     7.212673e-03     7.390252e-06
│    1551758188.5657811165    17     7.212673e-03     7.390252e-06
│    1551758223.2175629139    18     7.212639e-03     1.560990e-06
│    1551758257.7572259903    19     7.212608e-03     7.172192e-06
│    1551758290.823802948    20     7.212608e-03     5.154701e-08
│    1551758316.8794739246    21     7.212608e-03     5.154701e-08
│    1551758342.8025970459    22     7.212608e-03     5.154701e-08
│    1551758368.7285950184    23     7.212608e-03     5.154701e-08
│    1551758394.6267869473    24     7.212608e-03     5.154701e-08
│    1551758420.5754930973    25     7.212608e-03     5.154701e-08
│    1551758446.5254330635    26     7.212608e-03     5.154701e-08
│    1551758472.454957962    27     7.212608e-03     2.368872e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method FDDq, and printing time and q :
│
│    1.551758485966057e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1551758730.3310580254     0     6.521377e-02     5.500858e-02
│    1551758980.3469901085     1     1.388960e-02     2.473551e-02
│    1551759231.7757260799     2     7.817013e-03     2.440243e-03
│    1551759268.4188189507     3     7.817013e-03     2.440243e-03
│    1551759513.465708971     4     7.240150e-03     1.108362e-03
│    1551759545.0571510792     5     7.240150e-03     1.108362e-03
│    1551759790.3847711086     6     7.228540e-03     1.137964e-04
│    1551760027.7207050323     7     7.224389e-03     7.296822e-04
│    1551760264.1746320724     8     7.218334e-03     4.886249e-04
│    1551760502.5484640598     9     7.214718e-03     4.958740e-04
│    1551760740.1306209564    10     7.212923e-03     9.511154e-05
│    1551760976.5135378838    11     7.212629e-03     7.209988e-05
│    1551761210.6786019802    12     7.212608e-03     1.301734e-07
│    1551761236.66918993    13     7.212608e-03     1.301734e-07
│    1551761262.6539320946    14     7.212608e-03     1.301734e-07
│    1551761288.6356539726    15     7.212608e-03     1.301734e-07
│    1551761314.6406729221    16     7.212608e-03     1.301734e-07
│    1551761340.6277658939    17     7.212608e-03     1.301734e-07
│    1551761366.6241810322    18     7.212608e-03     1.301734e-07
│    1551761392.6452329159    19     7.212608e-03     1.301734e-07
│    1551761418.6504468918    20     7.212608e-03     1.147283e-08
└────────────────────────

┌────────────────────────
│ Optimizing using method ADDq, and printing time and q :
│
│    1.551761427202471e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1551761629.6983449459     0     6.521377e-02     5.500858e-02
│    1551761940.2861320972     1     1.367760e-02     2.434447e-02
│    1551762274.8026890755     2     7.791760e-03     2.402866e-03
│    1551762311.2957570553     3     7.791760e-03     2.402866e-03
│    1551762491.3001790047     4     7.238546e-03     1.034544e-03
│    1551762531.9407939911     5     7.238546e-03     1.034544e-03
│    1551763008.6284229755     6     7.228064e-03     1.213862e-04
│    1551763434.4954199791     7     7.223945e-03     7.629348e-04
│    1551763759.7934660912     8     7.217872e-03     4.117736e-04
│    1551763789.5299289227     9     7.217872e-03     4.117736e-04
│    1551764300.1290330887    10     7.214556e-03     4.384884e-04
│    1551764511.0417521    11     7.212904e-03     1.137828e-04
│    1551764789.6769599915    12     7.212623e-03     5.275763e-05
│    1551764815.7157950401    13     7.212623e-03     5.275763e-05
│    1551764841.7746610641    14     7.212623e-03     5.275763e-05
│    1551764867.7038359642    15     7.212623e-03     5.275763e-05
│    1551764893.7761321068    16     7.212623e-03     5.275763e-05
│    1551765233.697437048    17     7.212608e-03     8.153595e-07
│    1551765504.6984729767    18     7.212608e-03     2.052311e-08
│    1551765530.7620489597    19     7.212608e-03     1.092557e-11
└────────────────────────

┌────────────────────────
│ Optimizing using method CSDDq, and printing time and q :
│
│    1.551765539320366e9
│    time                   iteration      q(λ)           |Dq(λ)|
│    1551766004.9462089539     0     6.521377e-02     5.500858e-02
│    1551766644.1099970341     1     1.367760e-02     2.434447e-02
│    1551767213.5281050205     2     7.791760e-03     2.402866e-03
│    1551767250.170830965     3     7.791760e-03     2.402866e-03
│    1551767950.7097389698     4     7.238546e-03     1.034544e-03
│    1551767982.7712020874     5     7.238546e-03     1.034544e-03
│    1551768626.3271861076     6     7.228064e-03     1.213862e-04
│    1551769053.2136509418     7     7.223945e-03     7.629349e-04
│    1551769380.238038063     8     7.217872e-03     4.117735e-04
│    1551769410.0014419556     9     7.217872e-03     4.117735e-04
│    1551769956.7617139816    10     7.214556e-03     4.384883e-04
│    1551769985.4385509491    11     7.214556e-03     4.384883e-04
│    1551770404.5565180779    12     7.213443e-03     1.422530e-05
│    1551770434.0334420204    13     7.213443e-03     1.422530e-05
│    1551771309.8578929901    14     7.212996e-03     2.642094e-05
│    1551771816.9819290638    15     7.212656e-03     1.104680e-04
│    1551772265.233741045    16     7.212608e-03     1.281638e-08
│    1551772291.276966095    17     7.212608e-03     3.170371e-11
└────────────────────────

=#
