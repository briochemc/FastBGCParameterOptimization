# Set the options for the Newton optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = false, x_tol = 1e-3)

# Dictionary to hold the results
numFac = Dict()

# Using Cassette to plot benchmark of convergence
Cassette.@context BenchmarkData

# Add prehook for storing number of calls of factorize
Cassette.prehook(ctx::BenchmarkData, ::typeof(factorize), args...) = ctx.metadata.factorizations[] += 1
# Add posthook for storing values of objective and time at which it was computed
myreal(x::Dual) = DualNumbers.realpart(x)
myreal(x::Hyper) = HyperDualNumbers.realpart(x)
myreal(x::Complex) = real(x)
myreal(x::Float64) = x
function Cassette.posthook(ctx::BenchmarkData, output, ::typeof(q!), args...)
    push!(ctx.metadata.qvalues, myreal(output))
    push!(ctx.metadata.factorization_counter, ctx.metadata.factorizations[])
    push!(ctx.metadata.ftimer, time())
end


# Initiate BenchmarkData by creating an instance of the `tape`
struct ProfileCtx
    factorizations::Ref{Int64}
    factorization_counter::Vector{Int64}
    qvalues::Vector{Float64}
    ftimer::Vector{Float64}
end
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))

# make it into short code by listing the methods differently and
# interpolating them using the $ sign
#list_methods = [
#                (q!, Dq!, D2q!)
#               ]


init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, D2q!, λ₀, Newton(), opt))
numFac["D2q"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, CSDDq!, λ₀, Newton(), opt))
numFac["CSDDq"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, ADq!, AD2q!, λ₀, Newton(), opt))
numFac["AD2q"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, FDDq!, λ₀, Newton(), opt))
numFac["FDDq"] = tape.metadata

# Do it again because of the timer

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, D2q!, p2λ(3p₀), Newton(), opt))
numFac["D2q"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, CSDDq!, p2λ(3p₀), Newton(), opt))
numFac["CSDDq"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, ADq!, AD2q!, p2λ(3p₀), Newton(), opt))
numFac["AD2q"] = tape.metadata

init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
tape = BenchmarkData(metadata=ProfileCtx(Ref(0), [], [], []))
Cassette.@overdub(tape, optimize(q!, Dq!, FDDq!, p2λ(3p₀), Newton(), opt))
numFac["FDDq"] = tape.metadata



# # Now plot the results using `Plots.jl`
# using Plots
# 
# # Convergence vs time
# timer = tape.metadata.ftimer .- tape.metadata.ftimer[1]
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
