# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = false, x_tol = 1e-3)

# Using Cassette to plot benchmark of convergence
Cassette.@context BenchmarkData

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
function Cassette.prehook(ctx::BenchmarkData, output, ::typeof(q!), args...)
    ctx.metadata.q_counter[] += 1
    push!(ctx.metadata.qtimes_before, time())
    push!(ctx.metadata.q_calls, ctx.metadata.q_counter[])
end
# Add posthook for storing values of objective and time at which it was computed
myreal(x::Dual) = DualNumbers.realpart(x)
myreal(x::Hyper) = HyperDualNumbers.realpart(x)
myreal(x::Complex) = real(x)
myreal(x::Float64) = x
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

# Create the context type for storing the data
@default_kw struct ProfileCtx
    factorize_counter::Ref{Int64}          | Ref(0)
    f_counter::Ref{Int64}                  | Ref(0)
    q_counter::Ref{Int64}                  | Ref(0)
    factorize_calls::Vector{Int64}         | []
    f_calls::Vector{Int64}                 | []
    q_calls::Vector{Int64}                 | []
    factorize_type::Vector{String}         | []
    factorizetimes_before::Vector{Float64} | []
    factorizetimes_after::Vector{Float64}  | []
    ftimes_before::Vector{Float64}         | []
    ftimes_after::Vector{Float64}          | []
    qtimes_before::Vector{Float64}         | []
    qtimes_after::Vector{Float64}          | []
    qvalues::Vector{Float64}               | []
end

# Reformat the context type into a tuple of standard Julia types
read_profile(p::ProfileCtx) = Dict(
    "factorize_calls" => p.factorize_calls,
    "f_calls" => p.f_calls,
    "q_calls" => p.q_calls,
    "factorize_type" => p.factorize_type,
    "factorizetimes_before" => p.factorizetimes_before,
    "factorizetimes_after" => p.factorizetimes_after,
    "ftimes_before" => p.ftimes_before,
    "ftimes_after" => p.ftimes_after,
    "qtimes_before" => p.qtimes_before,
    "qtimes_after" => p.qtimes_after,
    "qvalues" => p.qvalues
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
