# Set the options for the NewtonTrustRegion optimizer
opt = Optim.Options(store_trace = true, show_trace = false, extended_trace = false, x_tol = 1e-3)

# run once to trigger compilation
init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
optimize(q!, Dq!, D2q!, λ₀, NewtonTrustRegion(), opt)

using Profile
Profile.clear()  # in case we have any previous profiling data
Profile.init(n = 10^7) # increase buffer size
init.x, init.p = 1x₀, 3p₀
J.fac, J.p = factorize(fJac(x₀, 3p₀)), 3p₀
@profile optimize(q!, Dq!, D2q!, λ₀, NewtonTrustRegion(), opt)

li, lidict = Profile.retrieve()
using JLD2

path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
jld_file = joinpath(path_to_package_root, "data/D2q_profile_data.jld2")

@save jld_file li lidict

# Other method to save profile output
r = Profile.retrieve();
profile_file = joinpath(path_to_package_root, "data/D2q_profile_data.bin")
profile_file_pointer = open(profile_file, "w")
using Serialization
serialize(profile_file_pointer, r)
close(profile_file_pointer)

# path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-1]...)
# jld_file = joinpath(path_to_package_root, "data/D2q_profile_data.jld2")
# @load jld_file
# ProfileView.view(li, lidict=lidict)

# Code to plot the profile
using ProfileView, Serialization
profile_file = "data/D2q_profile_data.bin"
f = open(profile_file)
r = deserialize(f);
ProfileView.view(r[1], lidict=r[2])

