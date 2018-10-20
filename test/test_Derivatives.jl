#=-
These tests should be run everytime to make sure the tests are passed on the functions defined by:
    - bgc functions
    - cost_functions
=#




# TODO Put them in Continuous Integration

# Test the Jacobian of f





@testset "Checks fJac(x, p) ≈ numJac(x, p)" begin
    C = build_big_compressor(wet3d, iwet, nwet, 2, 2)
    for i in 1:10, j in 1:10
        x = exp.(randn(length(x₀))) .* x₀
        p = Para(exp.(randn(length(vec(p₀))))) * p₀ # * overloaded to multiply element-wise Para types
        @test fJac(x, p) ≈ numJac(x -> f(x, p), x, C)
    end
end

# Test the gradient of q (not working yet)
"""
    ComplexStepGradient(q, λ::Vector{U})

Returns the gradient using the complex step method.
Only good for small sizes.
`q` is an application from Rⁿ to R in this case.
"""
function ComplexStepGradient(q, λ::Vector{U}) where U # q : Rⁿ -> R
    n = length(λ)
    out = zeros(U, n)
    h = 1e-50
    for i in 1:length(λ)
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[i] = imag.(q(imλ)) / h
    end
    return out
end
CSDq!(λ) = ComplexStepGradient(q!, λ)'

"""
    ComplexStepJacobian(Dq, λ::Vector{U})

Returns the Jacobian using the complex step method.
Only good for small sizes.
(Do not use for the state model `J` if using OCIM!)
`Dq` is an application from Rⁿ to Rⁿ in this case.
"""
function ComplexStepJacobian(Dq, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    h = 1e-50
    for i in 1:n
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[:, i] .= imag.(vec(Dq(imλ))) / h
    end
    return out
end
CSD2q!(λ) = ComplexStepJacobian(Dq!, λ)


@testset "Check Derivatives of cost functions" begin
    @testset "Check `Dq!`" begin
        @testset "Against Calculus.gradient (large rtol)" begin
            @test isapprox(gradient_q!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Against Calculus.jacobian (large rtol)" begin
            @test isapprox(jacobian_q!(λ₀), Dq!(λ₀), rtol=1e-3)
        end
        @testset "Against complex-step gradient (mine)" begin
            @test isapprox(CSDq!(λ₀), Dq!(λ₀))
        end
    end
    @testset "Check `D2q!`" begin
        @testset "Against Calculus.jacobian (large rtol)" begin
            @test isapprox(jacobian_Dq!(λ₀), D2q!(λ₀), rtol=1e-3)
        end
        @testset "Against complex-step jacobian (mine)" begin
            @test isapprox(CSD2q!(λ₀), D2q!(λ₀))
        end
    end
end

