# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = false, x_tol = 1e-3)

# Using Cassette to plot benchmark of convergence
Cassette.@context BenchmarkData

myreal(x::Dual) = DualNumbers.realpart(x)
myreal(x::Hyper) = HyperDualNumbers.realpart(x)
myreal(x::Complex) = real(x)
myreal(x::Float64) = x

functions_timed = [
    (:Dq,    :Dq!),
    (:FDq,   :FDq!),
    (:CSDq,  :CSDq!),
    (:ADq,   :ADq!),
    (:D2q,   :D2q!),
    (:FD2q,  :FD2q!),
    (:FDDq,  :FDDq!),
    (:CSDDq, :CSDDq!),
    (:ADDq,  :ADDq!),
    (:AD2q,  :AD2q!),
    (:factorize, :factorize),
    (:backsolve, :\)
]

for (fsym1, fsym2) in functions_timed
    @eval function Cassette.prehook(ctx::BenchmarkData, ::typeof($fsym2), args...)
        push!(ctx.metadata.$fsym1.tics, time())
    end
    @eval function Cassette.posthook(ctx::BenchmarkData, output, ::typeof($fsym2), args...)
        push!(ctx.metadata.$fsym1.tocs, time())
    end
end
function Cassette.posthook(ctx::BenchmarkData, output, ::typeof(q!), args...)
    push!(ctx.metadata.q.tocs, time())
    push!(ctx.metadata.qvalues, myreal(output))
end


#=

@default_kw struct ProfileCtx
    all_times::Timers         | Timers()
    q_values::Vector{Float64} | []
    f_norms::Vector{Float64}  | []
end

# Macro to generate Timers
@default_kw struct Timers

end


# Add prehooks for storing number of calls
function Cassette.prehook(ctx::BenchmarkData, ::typeof(factorize), args...)
    ctx.metadata.factorize_counter[] += 1
    push!(ctx.metadata.factorize_calls, ctx.metadata.factorize_counter[])
    push!(ctx.metadata.factorizetimes_before, time())
    push!(ctx.metadata.factorize_type, string(typeof(args[1])))
end
function Cassette.prehook(ctx::BenchmarkData, ::typeof(f), args...)
    ctx.metadata.f_counter[] += 1
    push!(ctx.metadata.ftimes_before, time())
    push!(ctx.metadata.f_calls, ctx.metadata.f_counter[])
end
function Cassette.prehook(ctx::BenchmarkData, ::typeof(q!), args...)
    ctx.metadata.q_counter[] += 1
    push!(ctx.metadata.qtimes_before, time())
    push!(ctx.metadata.q_calls, ctx.metadata.q_counter[])
end
# Add posthook for storing values of objective and time at which it was computed
function Cassette.posthook(ctx::BenchmarkData, output, ::typeof(q!), args...)
    push!(ctx.metadata.qvalues, myreal(output))
    push!(ctx.metadata.qtimes_after, time())
end
function Cassette.posthook(ctx::BenchmarkData, output, ::typeof(f), args...)
    push!(ctx.metadata.ftimes_after, time())
end
function Cassette.posthook(ctx::BenchmarkData, output, ::typeof(factorize), args...)
    push!(ctx.metadata.factorizetimes_after, time())
end
=#



@default_kw struct FunctionTimer
    tics::Vector{Float64} | []
    tocs::Vector{Float64} | []
end

# Create the context type for storing the data
@default_kw struct ProfileCtx
    qvalues::Vector{Float64} | []
    q::FunctionTimer | FunctionTimer()
    Dq::FunctionTimer | FunctionTimer()
    FDq::FunctionTimer | FunctionTimer()
    CSDq::FunctionTimer | FunctionTimer()
    ADq::FunctionTimer | FunctionTimer()
    D2q::FunctionTimer | FunctionTimer()
    FD2q::FunctionTimer | FunctionTimer()
    FDDq::FunctionTimer | FunctionTimer()
    CSDDq::FunctionTimer | FunctionTimer()
    ADDq::FunctionTimer | FunctionTimer()
    AD2q::FunctionTimer | FunctionTimer()
    factorize::FunctionTimer | FunctionTimer()
    backsolve::FunctionTimer | FunctionTimer()
end

# Reformat the context type into a tuple of standard Julia types
read_profile(p::ProfileCtx) = Dict(
    "qvalues" => p.qvalues,
    "q" => (p.q.tics, p.q.tocs),
    "Dq" => (p.Dq.tics, p.Dq.tocs),
    "FDq" => (p.FDq.tics, p.FDq.tocs),
    "CSDq" => (p.CSDq.tics, p.CSDq.tocs),
    "ADq" => (p.ADq.tics, p.ADq.tocs),
    "D2q" => (p.D2q.tics, p.D2q.tocs),
    "FD2q" => (p.FD2q.tics, p.FD2q.tocs),
    "FDDq" => (p.FDDq.tics, p.FDDq.tocs),
    "CSDDq" => (p.CSDDq.tics, p.CSDDq.tocs),
    "ADDq" => (p.ADDq.tics, p.ADDq.tocs),
    "AD2q" => (p.AD2q.tics, p.AD2q.tocs),
    "factorize" => (p.factorize.tics, p.factorize.tocs),
    "backsolve" => (p.backsolve.tics, p.backsolve.tocs)
)

# Dictionary to hold the results
# Load it if it exists, otherwise creat a new one
path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data/methods_benchmark_data.jld2")
isfile(jld_file) ? (@load jld_file methods_benchmark_data) : methods_benchmark_data = Dict()

# make it into short code by listing the methods differently and
# interpolating them using the $ sign
list_methods = [
    ( "D2q", :q!,  :Dq!,   :D2q!)
    ("AD2q", :q!, :ADq!,  :AD2q!)
    ("CSDq", :q!,  :Dq!, :CSDDq!)
    ("FDDq", :q!,  :Dq!,  :FDDq!)
]

for (method_name, q, Dq, D2q) in list_methods
    println("\n\n\n------------------------\n\n\n")
    println("Running " * method_name * "with Cassette:")
    println("\n\n\n------------------------\n\n\n")
    init.x, init.p = 1x₀, 3p₀
    J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
    local tape
    tape = BenchmarkData(metadata=ProfileCtx())
    eval( :( Cassette.@overdub($tape, optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)) ) )
    methods_benchmark_data[method_name * "_1strun"] = read_profile(tape.metadata)
    @save jld_file methods_benchmark_data
end
# Do it again because of the timer
for (method_name, q, Dq, D2q) in list_methods
    println("\n\n\n------------------------\n\n\n")
    println("Running " * method_name * "with Cassette:")
    println("\n\n\n------------------------\n\n\n")
    init.x, init.p = 1x₀, 3p₀
    J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
    tape = BenchmarkData(metadata=ProfileCtx())
    eval( :( Cassette.@overdub($tape, optimize($q, $Dq, $D2q, $λ₀, NewtonTrustRegion(), $opt)) ) )
    methods_benchmark_data[method_name * "_2ndrun"] = read_profile(tape.metadata)
    @save jld_file methods_benchmark_data
end






# # Now plot the results using `Plots.jl`
# using Plots
#
# # Convergence vs time
# timer = tape.metadata.qtimes .- tape.metadata.qtimes[1]
# y = tape.metadata.fvalues
# p1 = plot(timer, y, yaxis = :log)
# xlabel!("computing time (ms)")
# ylabel!("f")
# legend!("")
#
# # Convergence vs number of calls
# counter = tape.metadata.fcounter
# p2 = plot(counter, y, yaxis = :log)
# xlabel!("number of f calls")
# ylabel!("f")
#
# # Combine subplotds into single figure
# plot(p1, p2, layout = (2, 1))
