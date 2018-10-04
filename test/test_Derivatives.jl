#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#

# TODO Put them in Continuous Integration

# Test the Jacobian of f

@testset "Checks fJac(x, p) ≈ numjac(x, p)" begin
    C = build_big_compressor(wet3d, iwet, nwet, 2, 2)
    for i in 1:10, j in 1:10
        x = exp.(randn(length(x₀))) .* x₀
        p = Para(exp.(randn(length(vec(p₀))))) * p₀ # * overloaded to multiply element-wise Para types
        @test fJac(x, p) ≈ numJac(x -> f(x, p), x, C)
    end
end
@test numjac(x, p) = fJac(x, p)
foo = fJac(x₀, p₀)

