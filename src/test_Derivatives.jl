#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#

# TODO Put them in Continuous Integration

# Test the Jacobian of f
using Test
@testset "Checks fJac(x, p) ≈ numjac(x, p) $i $j" for i in 1:10, j in 1:10
    C = build_big_compressor(wet3d, iwet, nwet, 2, 2)
    x = exp.(rand(length(x₀))) .* x₀
    p = exp.(rand(length(p₀))) .* p₀
    @test fJac(x, p) ≈ numjac(x -> f(x, p), x, C)
end
@test numjac(x, p) = fJac(x, p)
foo = fJac(x₀, p₀)

