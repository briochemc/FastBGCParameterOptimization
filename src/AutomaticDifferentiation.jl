"""
    ComplexStepGradient(f̂, λ::Vector{U})

Returns the gradient using the complex step method.
Only good for small sizes.
`f̂` is an application from Rⁿ to R in this case.
"""
function ComplexStepGradient(f̂, λ::Vector{U}) where U # f̂ : Rⁿ -> R
    n = length(λ)
    out = zeros(U, n)
    h = 1e-50
    for i in 1:length(λ)
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[i] = imag.(f̂(imλ)) / h
    end
    return out
end

"""
    ComplexStepJacobian(∇f̂, λ::Vector{U})

Returns the Jacobian using the complex step method.
Only good for small sizes.
(Do not use for the state model `J` if using OCIM!)
`∇f̂` is an application from Rⁿ to Rⁿ in this case.
"""
function ComplexStepJacobian(∇f̂, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    h = 1e-50
    for i in 1:n
        imλ = convert(Vector{Complex{U}}, λ)
        imλ[i] += h * im
        out[:, i] .= imag.(vec(∇f̂(imλ))) / h
    end
    return out
end

"""
    DualNumbersGradient(f̂, λ::Vector{U})

Returns the gradient using DualNumbers.
Only good for small sizes.
`f̂` is an application from Rⁿ to R in this case.
"""
function DualNumbersGradient(f̂, λ::Vector{U}) where U # f̂ : Rⁿ -> R
    n = length(λ)
    out = zeros(U, n)
    for i in 1:length(λ)
        ελ = convert(Vector{Dual{U}}, λ)
        ελ[i] += ε
        out[i] = dualpart(f̂(ελ))
    end
    return out
end

"""
    DualNumbersJacobian(∇f̂, λ::Vector{U})

Returns the Jacobian using DualNumbers.
Only good for small sizes.
(Do not use for the state model `J` if using OCIM!)
`∇f̂` is an application from Rⁿ to Rⁿ in this case.
"""
function DualNumbersJacobian(∇f̂, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    for i in 1:n
        ελ = convert(Vector{Dual{U}}, λ)
        ελ[i] += ε
        out[:, i] .= dualpart.(vec(∇f̂(ελ)))
    end
    return out
end

"""
    HyperDualNumbersHessian(f̂, λ::Vector{U})

Returns the Hessian using HyperDualNumbers.
Only good for small sizes.
`f̂` is an application from Rⁿ to R in this case.
"""
function HyperDualNumbersHessian(f̂, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    for i in 1:n, j in 1:n
        hλ = convert(Vector{Hyper{U}}, λ)
        hλ[i] += ε₁
        hλ[j] += ε₂
        out[i, j] = ε₁ε₂part(f̂(hλ))
    end
    return out
end

"""
    HyperDualNumbersSymmetricHessian(f̂, λ::Vector{U})

Returns the Hessian using HyperDualNumbers.
Does not compute all the terms (because output should be symmetric)
Only good for small sizes.
`f̂` is an application from Rⁿ to R in this case.
"""
function HyperDualNumbersSymmetricHessian(f̂, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    for i in 1:n, j in i:n
        hλ = convert(Vector{Hyper{U}}, λ)
        hλ[i] += ε₁
        hλ[j] += ε₂
        out[i, j] = ε₁ε₂part(f̂(hλ))
        i ≠ j ? out[j, i] = out[i, j] : nothing
    end
    return out
end


function update_HD_buffer!(f̂, buffer, last_λ, λ::Vector{U}) where U<:Float64
    if λ != last_λ
        copyto!(last_λ, λ)
        n = length(λ)
        for i in 1:n
            hλ = convert(Vector{Hyper{U}}, λ)
            hλ[i] += ε₁
            hλ[1] += ε₂
            buffer[i] = f̂(hλ)
        end
    end
end

function HyperDualNumbersBufferedHessian!(f̂, buffer, last_λ, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    update_HD_buffer!(f̂, buffer, last_λ, λ::Vector{U})
    out[1, :] .= ε₁ε₂part.(buffer)
    out[:, 1] .= out[1, :]
    for i in 2:n, j in i:n
        hλ = convert(Vector{Hyper{U}}, λ)
        hλ[i] += ε₁
        hλ[j] += ε₂
        out[i, j] = ε₁ε₂part(f̂(hλ))
        i ≠ j ? out[j, i] = out[i, j] : nothing
    end
    return out
end

function HyperDualNumbersBufferedgradient!(f̂, buffer, last_λ, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    update_HD_buffer!(f̂, buffer, last_λ, λ::Vector{U})
    return ε₁part.(buffer)
end

Base.conj(x::Hyper) = Hyper(conj(real(x)), conj(ε₁part(x)), conj(ε₂part(x)), conj(ε₁ε₂part(x)))
