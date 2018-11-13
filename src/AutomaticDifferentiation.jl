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

"""
    DualNumbersGradient(q, λ::Vector{U})

Returns the gradient using DualNumbers.
Only good for small sizes.
`q` is an application from Rⁿ to R in this case.
"""
function DualNumbersGradient(q, λ::Vector{U}) where U # q : Rⁿ -> R
    n = length(λ)
    out = zeros(U, n)
    for i in 1:length(λ)
        ελ = convert(Vector{Dual{U}}, λ)
        ελ[i] += ε
        out[i] = dualpart.(q(ελ))
    end
    return out
end

"""
    DualNumbersJacobian(Dq, λ::Vector{U})

Returns the Jacobian using DualNumbers.
Only good for small sizes.
(Do not use for the state model `J` if using OCIM!)
`Dq` is an application from Rⁿ to Rⁿ in this case.
"""
function DualNumbersJacobian(Dq, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    for i in 1:n
        ελ = convert(Vector{Dual{U}}, λ)
        ελ[i] += ε
        out[:, i] .= dualpart.(vec(Dq(ελ)))
    end
    return out
end

"""
    DualNumbersHessian(q, λ::Vector{U})

Returns the Hessian using HyperDualNumbers.
Only good for small sizes.
`q` is an application from Rⁿ to R in this case.
"""
function DualNumbersHessian(Dq, λ::Vector{U}) where U<:Float64
    n = length(λ)
    out = zeros(U, n, n)
    for i in 1:n, j in 1:n
        hλ = convert(Vector{Hyper{U}}, λ)
        hλ[i] += ε₁ ## unfinished business here
        hλ[j] += ε₂ ## unfinished business here
        out[i, j] .= eps1eps2.(q(hλ))
    end
    return out
end