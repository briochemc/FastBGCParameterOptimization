#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#




# TODO Put them in Continuous Integration

# Test the Jacobian of f





#@testset "Checks fJac(x, p) ≈ numJac(x, p)" begin
#    C = build_big_compressor(wet3d, Circulation.iwet, nwet, 2, 2)
#    for i in 1:10
#        x = exp.(randn(length(x₀))) .* x₀
#        p = Para(exp.(randn(length(vec(p₀))))) * p₀ # * overloaded to multiply element-wise Para types
#        @test ∇ₓF(x, p) ≈ numJac(x -> f(x, p), x, C)
#    end
#end

Hessian_methods = [
    :OF1_∇²f̂!,
    :F1_∇²f̂!,
    :F0_∇²f̂!,
    :DUAL_∇²f̂!,
    :CSD_∇²f̂!,
    :HYPER_∇²f̂!
]

FD_Hessian_methods = [
    :FD1_∇²f̂!,
    :FD2_∇²f̂!
]

gradient_methods = [
    :∇f̂!,
    :OF1_∇f̂!,
    :F1_∇f̂!,
    :CSD_∇f̂!,
    :HYPER_∇f̂!
]

FD_gradient_methods = [
    :FD2_∇f̂!
]


λ = exp(randn()) .* λ₀

@testset "Check Derivatives of objective" begin
    @testset "Accurate methods" begin
        @testset "gradient" begin
            #end
            @testset "$m1 VS $m2" for m1 in gradient_methods, m2 in gradient_methods
                @eval @test isapprox($m1(λ), $m2(λ))
            end
        end
        @testset "Hessian" begin
            @testset "$m1 VS $m2" for m1 in Hessian_methods, m2 in Hessian_methods
                @eval @test isapprox($m1(λ), $m2(λ))
            end
        end
    end
    @testset "Finite-difference methods" begin
        @testset "gradient" begin
            @testset "Analytical VS Calculus.gradient (large rtol)" begin
                @test isapprox(gradient_f̂!(λ), ∇f̂!(λ), rtol=1e-3)
            end
            @testset "$m VS OF1" for m in FD_gradient_methods
                @eval @test isapprox($m(λ), OF1_∇f̂!(λ), rtol=1e-3)
            end
        end
        @testset "Hessian" begin
            @testset "$m VS OF1" for m in FD_Hessian_methods
                @eval @test isapprox($m(λ), OF1_∇²f̂!(λ), rtol=1e-3)
            end
        end
    end
end

