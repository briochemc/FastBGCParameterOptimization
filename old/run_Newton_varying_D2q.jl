# Set the options for the Newton optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = true, x_tol = 1e-3)

# Dictionary to hold the results
results = Dict()

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
@btime results["D2q"] = optimize(q!, Dq!, D2q!, λ₀, Newton(), opt)

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
@btime results["CSDDq"] = optimize(q!, Dq!, CSDDq!, λ₀, Newton(), opt)

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
@btime results["AD2q"] = optimize(q!, ADq!, AD2q!, λ₀, Newton(), opt)

init.x, init.p = x₀, p₀
J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
@btime results["FDDq"] = optimize(q!, Dq!, FDDq!, λ₀, Newton(), opt)

#init.x, init.p = x₀, p₀
#J.fac, J.p = factorize(fJac(x₀, p₀)), p₀
#@btime results["FD2q"] = optimize(q!, FDq!, FD2q!, λ₀, Newton(), opt)


