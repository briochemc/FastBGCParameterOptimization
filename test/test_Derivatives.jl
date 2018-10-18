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

function ComplexStepGradient(f, x::Vector{U}) where U # f : Rⁿ -> R
    n = length(x)
    out = zeros(U, n)
    h = 1e-50
    for i in 1:length(x)
        imx = convert(Vector{Complex{U}}, x)
        imx[i] += h * im
        out[i] = imag.(f(imx)) / h
    end
    return out
end


function ComplexStepJacobian(f, x::Vector{U}) where U # f : Rⁿ -> Rⁿ
    n = length(x)
    out = zeros(U, n, n)
    h = 1e-50
    for i in 1:n
        imx = convert(Vector{Complex{U}}, x)
        imx[i] += h * im
        out[:, i] .= imag.(vec(f(imx))) / h
    end
    return out
end
#
#@testset "Check that Dq!(λ) ≈ CSD(q!, λ)" begin
#    for i in 1:10
#        λ = randn(npopt)
#        @test vec(Dq!(λ)) ≈ CSD(q!, λ)
#    end
#end
#
#@testset "Check that D2q!(λ) ≈ Calculus.Gradient(Dq!, λ)" begin
#    for i in 1:10
#        λ = randn(npopt)
#        @test vec(Dq!(λ)) ≈ Calculus.gradient(q!, λ) rtol = 1e-3
#    end
#end
