# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = false, x_tol = 1e-3)

# run once to trigger compilation
init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
optimize(q!, Dq!, D2q!, λ₀, NewtonTrustRegion(), opt)

using Profile
Profile.clear()  # in case we have any previous profiling data
init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
@profile optimize(q!, Dq!, D2q!, λ₀, NewtonTrustRegion(), opt)

using ProfileView

ProfileView.view()
