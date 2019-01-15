# Set the options for the Newton optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = true, x_tol = 1e-3)

# Dictionary to hold the results
numFac = Dict()

# Using Cassette to count the number of factorizations
# Context for counting factorizations
Cassette.@context CountFactorizations
Cassette.prehook(ctx::CountFactorizations, ::typeof(factorize), args...) = ctx.metadata[] += 1

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
counter = CountFactorizations(metadata=Ref(0))
Cassette.@overdub(counter, optimize(q!, Dq!, D2q!, λ₀, Newton(), opt))
numFac["D2q"] = counter.metadata.x

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
counter = CountFactorizations(metadata=Ref(0))
Cassette.@overdub(counter, optimize(q!, Dq!, CSDDq!, λ₀, Newton(), opt))
numFac["CSDDq"] = counter.metadata.x

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
counter = CountFactorizations(metadata=Ref(0))
Cassette.@overdub(counter, optimize(q!, ADq!, AD2q!, λ₀, Newton(), opt))
numFac["AD2q"] = counter.metadata.x

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
counter = CountFactorizations(metadata=Ref(0))
Cassette.@overdub(counter, optimize(q!, Dq!, FDDq!, λ₀, Newton(), opt))
numFac["FDDq"] = counter.metadata.x
