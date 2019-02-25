using Plots
plotlyjs()



TimerOutputs = [
  # CSDDq  # ADDq    # FDDq      # D2q
    3116   2489      1148         77.3       # D2q
     242   243        231        229         # Dq
]


p = plot(transpose(TimerOutputs), seriestype=:bar,
    labels=["D2q" "Dq"], legend=:top,
    xticks=(1:4, ["Complex step", "Dual numbers", "Finite difference", "This paper"]),
    xlabel="Method",
    ylabel="Time (seconds)",
    yminorticks=collect(0:60:3600),
    yticks=collect(0:300:3600)
)

savefig(p, "fig/TimerOutputs_data_katana.pdf")


#=

DATA FROM KATANA Feb 25 2019

 ───────────────────────────────────────────────────────────────────────
                                Time                   Allocations
                        ──────────────────────   ───────────────────────
    Tot / % measured:        2732s / 100%             858GiB / 100%

 Section        ncalls     time   %tot     avg     alloc   %tot      avg
 ───────────────────────────────────────────────────────────────────────
 Trust Region        1    2732s   100%   2732s    858GiB  100%    858GiB
   D2q (dual)        5    2489s  91.1%    498s    786GiB  91.6%   157GiB
   Dq                6     243s  8.90%   40.5s   71.9GiB  8.38%  12.0GiB
   q                 6   17.4ms  0.00%  2.90ms   45.9MiB  0.01%  7.64MiB
 ───────────────────────────────────────────────────────────────────────


                                   Time                   Allocations
                           ──────────────────────   ───────────────────────
     Tot / % measured:          3359s / 100%             870GiB / 100%

 Section           ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────
 Trust Region           1    3359s   100%   3359s    870GiB  100%    870GiB
   D2q (complex)        5    3116s  92.8%    623s    799GiB  91.8%   160GiB
   Dq                   6     242s  7.21%   40.4s   71.7GiB  8.24%  12.0GiB
   q                    6   16.9ms  0.00%  2.82ms   45.9MiB  0.01%  7.64MiB
 ──────────────────────────────────────────────────────────────────────────


                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:            1379s / 100%             385GiB / 100%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Trust Region               1    1379s   100%   1379s    385GiB  100%    385GiB
   D2q (finite-diff)        5    1148s  83.2%    230s    316GiB  82.1%  63.2GiB
   Dq                       6     231s  16.8%   38.5s   69.0GiB  17.9%  11.5GiB
   q                        6   16.6ms  0.00%  2.76ms   45.9MiB  0.01%  7.64MiB
 ──────────────────────────────────────────────────────────────────────────────

 ───────────────────────────────────────────────────────────────────────
                                Time                   Allocations
                        ──────────────────────   ───────────────────────
    Tot / % measured:         306s / 100%            90.7GiB / 100%

 Section        ncalls     time   %tot     avg     alloc   %tot      avg
 ───────────────────────────────────────────────────────────────────────
 Trust Region        1     306s   100%    306s   90.7GiB  100%   90.7GiB
   Dq                6     229s  74.8%   38.2s   69.2GiB  76.3%  11.5GiB
   D2q (mine)        5    77.3s  25.2%   15.5s   21.5GiB  23.7%  4.29GiB
   q                 6   16.6ms  0.01%  2.77ms   45.9MiB  0.05%  7.64MiB
 ───────────────────────────────────────────────────────────────────────


=#
