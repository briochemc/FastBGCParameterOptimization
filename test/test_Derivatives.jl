#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#




# TODO Put them in Continuous Integration

# Test the Jacobian of f





@testset "Checks fJac(x, p) ≈ numJac(x, p)" begin
    C = build_big_compressor(wet3d, Circulation.iwet, nwet, 2, 2)
    for i in 1:10
        x = exp.(randn(length(x₀))) .* x₀
        p = Para(exp.(randn(length(vec(p₀))))) * p₀ # * overloaded to multiply element-wise Para types
        @test fJac(x, p) ≈ numJac(x -> f(x, p), x, C)
    end
end

# Test the gradient of q (not working yet)



@testset "Check Derivatives of objective" begin
    @testset "gradient" begin
        @testset "Analytical VS Calculus.gradient (large rtol)" begin
            @test isapprox(gradient_q!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Analytical VS Calculus.jacobian (large rtol)" begin
            @test isapprox(FDq!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Analytical VS complex-step gradient (mine)" begin
            @test isapprox(CSDq!(λ₀), Dq!(λ₀))
        end
        @testset "Analytical VS DualNumbers gradient (mine)" begin
            @test isapprox(ADq!(λ₀), Dq!(λ₀))
        end
        @testset "Analytical VS HyperDualNumbers gradient (mine)" begin
            @test isapprox(HSDq!(λ₀), Dq!(λ₀))
        end
    end
    @testset "Check `D2q!`" begin
        @testset "FLASH VS Calculus.jacobian (large rtol)" begin
            @test isapprox(FDDq!(λ₀), D2q!(λ₀), rtol=1e-3)
        end
        @testset "FLASH VS Calculus.hessian (large rtol)" begin
            @test isapprox(FD2q!(λ₀), D2q!(λ₀), rtol=1e-3)
        end
        @testset "FLASH VS CSD" begin
            @test isapprox(CSDDq!(λ₀), D2q!(λ₀))
        end
        @testset "FLASH VS DUAL" begin
            @test isapprox(ADDq!(λ₀), D2q!(λ₀))
        end
        @testset "FLASH VS HYPER" begin
            @test isapprox(AD2q!(λ₀), D2q!(λ₀))
        end
        @testset "FLASH VS HYPERSMART" begin
            @test isapprox(HSD2q!(λ₀), D2q!(λ₀))
        end
    end
end

