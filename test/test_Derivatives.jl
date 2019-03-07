#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#




# TODO Put them in Continuous Integration

# Test the Jacobian of f





@testset "Checks fJac(x, p) ≈ numJac(x, p)" begin
    C = build_big_compressor(wet3d, SixBoxModel.iwet, nwet, 2, 2)
    for i in 1:10, j in 1:10
        x = exp.(randn(length(x₀))) .* x₀
        p = Para(exp.(randn(length(vec(p₀))))) * p₀ # * overloaded to multiply element-wise Para types
        @test fJac(x, p) ≈ numJac(x -> f(x, p), x, C)
    end
end

# Test the gradient of q (not working yet)



@testset "Check Derivatives of cost functions" begin
    @testset "Check `Dq!`" begin
        @testset "Against Calculus.gradient (large rtol)" begin
            @test isapprox(gradient_q!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Against Calculus.jacobian (large rtol)" begin
            @test isapprox(FDq!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Against complex-step gradient (mine)" begin
            @test isapprox(CSDq!(λ₀), Dq!(λ₀))
        end
        @testset "Against DualNumbers gradient (mine)" begin
            @test isapprox(ADq!(λ₀), Dq!(λ₀))
        end
    end
    @testset "Check `D2q!`" begin
        @testset "Against Calculus.jacobian (large rtol)" begin
            @test isapprox(FDDq!(λ₀), D2q!(λ₀), rtol=1e-3)
        end
        @testset "Against Calculus.hessian (large rtol)" begin
            @test isapprox(FD2q!(λ₀), D2q!(λ₀), rtol=1e-3)
        end
        @testset "Against complex-step jacobian (mine)" begin
            @test isapprox(CSDDq!(λ₀), D2q!(λ₀))
        end
        @testset "Against DualNumbers jacobian (mine)" begin
            @test isapprox(ADDq!(λ₀), D2q!(λ₀))
        end
        @testset "Against HyperDualNumbers Hessian (mine)" begin
            @test isapprox(AD2q!(λ₀), D2q!(λ₀))
        end
    end
end

