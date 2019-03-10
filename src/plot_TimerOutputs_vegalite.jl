# Load plotting package
using VegaLite, DataFrames # for stacked bars






# using LaTeXStrings # for LaTeX labels
#using ColorBrewer # for better colors


# Load data
using JLD2
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
str_out = "_default"
jld_file = joinpath(path_to_package_root, "data", "TimerOutputs_data" * str_out * "_katana" * ".jld2")
@load jld_file timers

translate_for_legend = Dict("q"=>"objective", "Dq"=>"gradient", "D2q"=>"Hessian")

# Reshape data into a DataFrame
function timer_to_DataFrame(timers::Dict)
    df = DataFrame(
        method = Array{String}(undef, 0),
        fgh = Array{String}(undef, 0),
        time = Array{Float64}(undef, 0),
        allocs = Array{Float64}(undef, 0),
        ncalls = Array{Int64}(undef, 0)
    )
    for k in keys(timers)
        if ~occursin("compiling", string(k))
            t = timers[k]
            for k2 in ["q", "Dq", "D2q"]
                push!(
                    df,
                    Dict(
                        :method => string(k),
                        :fgh => translate_for_legend[k2],
                        :time => t[k2]["time"] * 1e-9,
                        :allocs => t[k2]["allocs"] * 2^-30,
                        :ncalls => t[k2]["ncalls"]
                    )
                )
            end
        end
    end
    return df
end
df = timer_to_DataFrame(timers)

ptime = df |> @vlplot(
    :bar,
    x={"sum(time)", title="Computation time (seconds)"},
    y={:method, sort=["FLASH", "Dual", "FiniteDiff", "Complex", "HyperDual"]},
    color={:fgh, title="", scale={scheme="set2"}, sort=["q", "dq", "d²q"]},
    order={field=:method, sort=["q", "dq", "d²q"]}
)

# print to PDF
pdf_file_time = joinpath(path_to_package_root, "fig", "TimerOutputs_time.pdf")
save(pdf_file_time, ptime)


pmem = df |> @vlplot(
    :bar,
    x={"sum(allocs)", title="Memory allocations (GiB)"},
    y=:method,
    color={:fgh, title=""}
)

# print to PDF
pdf_file_mem = joinpath(path_to_package_root, "fig", "TimerOutputs_mem.pdf")
save(pdf_file_mem, pmem)

pncalls = df |> @vlplot(
    :bar,
    x={"sum(ncalls)", title="Number of calls"},
    y=:method,
    color={:fgh, title=""}
)

# print to PDF
pdf_file_ncalls = joinpath(path_to_package_root, "fig", "TimerOutputs_ncalls.pdf")
save(pdf_file_ncalls, pncalls)

#savefig(p, "fig/TimerOutputs_data_test.pdf")
#
#display(p)

## Data from printed output copied below
#TimerOutputs = [
#  # CSDDq  # ADDq       # FDDq         # D2q
#     229.0    223.0       214.0         223.0        # Dq
#    2978.0   1551.0      1065.0         81.3       # D2q
#]
#
#p = groupedbar(
#    color=ColorBrewer.palette("Set2", 3)[1:2]',
#    TimerOutputs',
#    bar_position=:stack,
#    labels=["Gradient", "Hessian"],
#    legend=:top,
#    xticks=(1:4, ["Complex step", "Dual numbers", "Finite difference", "This paper"]),
#    xlabel="Method",
#    ylabel="Time (seconds)",
#    yminorticks=collect(0:60:3600),
#    yticks=collect(0:300:3600)
#)
#
#display(p)
#
#savefig(p, "fig/TimerOutputs_data_katana.pdf")
#
#display(p)
#=

DATA FROM KATANA Feb 26 2019

 ───────────────────────────────────────────────────────────────────────
                                Time                   Allocations
                        ──────────────────────   ───────────────────────
    Tot / % measured:        1773s / 100%             540GiB / 100%

 Section        ncalls     time   %tot     avg     alloc   %tot      avg
 ───────────────────────────────────────────────────────────────────────
 Trust Region        1    1773s   100%   1773s    540GiB  100%    540GiB
   D2q (dual)        5    1551s  87.4%    310s    471GiB  87.3%  94.2GiB
   Dq                6     223s  12.6%   37.1s   68.7GiB  12.7%  11.5GiB
   q                 6   15.6ms  0.00%  2.61ms   45.8MiB  0.01%  7.64MiB
 ─────────────────────────────────────────────────────────────────────── ──────────────────────────────────────────────────────────────────────────
                                   Time                   Allocations
                           ──────────────────────   ───────────────────────
     Tot / % measured:          3207s / 100%             868GiB / 100%

 Section           ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────
 Trust Region           1    3207s   100%   3207s    868GiB  100%    868GiB
   D2q (complex)        5    2978s  92.9%    596s    797GiB  91.8%   159GiB
   Dq                   6     229s  7.14%   38.2s   71.2GiB  8.20%  11.9GiB
   q                    6   15.3ms  0.00%  2.55ms   45.8MiB  0.01%  7.64MiB
 ────────────────────────────────────────────────────────────────────────── ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:            1278s / 100%             384GiB / 100%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Trust Region               1    1278s   100%   1278s    384GiB  100%    384GiB
   D2q (finite-diff)        5    1065s  83.3%    213s    316GiB  82.1%  63.1GiB
   Dq                       6     214s  16.7%   35.6s   68.6GiB  17.8%  11.4GiB
   q                        6   13.9ms  0.00%  2.32ms   45.8MiB  0.01%  7.64MiB
 ──────────────────────────────────────────────────────────────────────────────

 ───────────────────────────────────────────────────────────────────────
                                Time                   Allocations
                        ──────────────────────   ───────────────────────
    Tot / % measured:         304s / 100%            90.2GiB / 100%

 Section        ncalls     time   %tot     avg     alloc   %tot      avg
 ───────────────────────────────────────────────────────────────────────
 Trust Region        1     304s   100%    304s   90.2GiB  100%   90.2GiB
   Dq                6     223s  73.3%   37.2s   68.7GiB  76.2%  11.5GiB
   D2q (mine)        5    81.3s  26.7%   16.3s   21.5GiB  23.8%  4.29GiB
   q                 6   15.6ms  0.01%  2.60ms   45.8MiB  0.05%  7.64MiB
 ───────────────────────────────────────────────────────────────────────


=#
