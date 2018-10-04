## This script ensures that there is a local copy of the unregistered packages.
using Pkg
Pkg.update()
Pkg.add(PackageSpec(url="https://github.com/rafaqz/Defaults.jl"))
Pkg.add(PackageSpec(url="https://github.com/briochemc/TransportMatrixTools.jl"))
