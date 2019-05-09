
#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
const S1 = buildPFD(ones(nb), DIV, Iabove)
const Sz = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S1 + w′ * Sz
end
T_all = (T_DIP, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(DIP, p)
    τg, DIPgeo = p.τg, p.DIPgeo
    return (DIPgeo .- DIP) ./ τg
end
# Uptake of phosphate (DIP)
relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    return Umax * relu(DIP) ./ (relu(DIP) .+ ku) .* (z .≤ z₀)
end
# Remineralization of particulate organic phosphorus (POP)
function remineralization(POP, p)
    κ = p.κ
    return κ * POP
end
# Add them up into sms functions (Sources Minus Sinks)
sms_DIP(DIP, POP, p) = geores(DIP, p) .- uptake(DIP, p) .+ remineralization(POP, p)
sms_POP(DIP, POP, p) = uptake(DIP, p) .- remineralization(POP, p)
sms_all = (sms_DIP, sms_POP) # bundles all the source-sink functions in a tuple

#===========================================
Let AIBECS generate 
the state function, F, and
the sparse Jacobian, ∇ₓF
===========================================#
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb) # generates the state function (and its Jacobian!)

#===========================================
Analytical derivatives
Required for CSD, DUAL, and HYPER methods
===========================================#

# Geological restoring
function geores(DIP, p)
    τg = p.τg, DIPgep = p.DIPgeo
    return (DIPgeo .- DIP) ./ τg
end
function georesJac(p)
    τg = p.τg
    return -I / τg
end
function ∂geores_∂DIPgeo(DIP, p)
    τg = p.τg, DIPgep = p.DIPgeo
    return 1 / τg
end

# Uptake
relu(x) = (x .≥ 0) .* x
drelu(x) = (x .≥ 0) .* 1.0
# Michaelis-Menten
mm(x, μ, k)     =  μ * x ./ (x .+ k)
∂mm_∂x(x, μ, k) =  μ * k ./ (x .+ k).^2
∂mm_∂μ(x, μ, k) =      x ./ (x .+ k)
∂mm_∂k(x, μ, k) = -μ * x ./ (x .+ k).^2
function uptake(DIP, p)
    umax, ku = p.umax, p.ku
    return d₀(maskEup) * mm(relu(DIP), umax, ku)
end
# Uptake derivatives
function uptakeJac(DIP, p)
    umax, ku = p.umax, p.ku
    return d₀(maskEup .* ∂mm_∂x(relu(DIP), umax, ku) .* drelu(DIP))
end
function ∂uptake_∂umax(DIP, p)
    umax, ku = p.umax, p.ku
    return maskEup .* ∂mm_∂μ(relu(DIP), umax, ku)
end
function ∂uptake_∂ku(DIP, p)
    umax, ku = p.umax, p.ku
    return maskEup .* ∂mm_∂k(relu(DIP), umax, ku)
end

# Remineralization
function remineralization(POP, p)
    κ = p.κ
    return κ * POP
end
# Remineralization derivatives
function remineralizationJac(POP, p)
    κ = p.κ
    return κ * I
end
function ∂remineralization_∂κ(POP, p)
    κ = p.κ
    return POP
end

# Indices for DIP and POP (needed for Rate of change F(x,p))
const iDIP = 1:nb
const iPOP = iDIP .+ nb
function unpackx(x)
    DIP = x[iDIP]
    POP = x[iPOP]
    return DIP, POP
end
# add method to deal with arrays for numJac
function unpackx(x::Array{<:Number,2})
    DIP = x[iDIP,:]
    POP = x[iPOP,:]
    return DIP, POP
end

# PFD transport (needed for Rate of change F(x,p))
const S1 = buildPFD(ones(nb), DIV, Iabove)
const Sz = buildPFD(ztop, DIV, Iabove)
function S(p::Para)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S1 + w′ * Sz
end

# Rate of change F(x,p)
function F(x, p)
    DIP, POP = unpackx(x)
    u = uptake(DIP, p)
    r = remineralization(POP, p)
    return [    -T * DIP - u + r + geores(DIP, p) ;
             -S(p) * POP + u - r                  ]
end


# Jacobian of f with respect to x
function ∇ₓF(x, p)
    DIP, POP = unpackx(x)
    uJac = uptakeJac(DIP, p)
    rJac = remineralizationJac(POP, p)
    foo = [ -T   - uJac + georesJac(p)    rJac        ;
                  +uJac                  -rJac - S(p) ]
    dropzeros!(foo)
    return foo
end

using Flatten
"""
    ∇ₚF(x, p)

Evaluates the jacobian of `f` with respect to `p`.
Concatenates `∇ₚF(x, p::Para, s::Symbol)` for all optimizable parameter symbols `s`.
"""
∇ₚF(x, p) = hcat((∇ₚF(x, p, s) for s in fieldnameflatten(p))...)

"""
    ∇ₚF(x, p, s::Symbol)

Evaluates the derivative of `F` with respect to `p.s` where `s` is the name of the parameter (of type `Symbol`).

You should fill this function with all the first derivatives of f with respoect to each parameter.
Called without the symbol `s`, this function will loop through all the optimizable parameters and create the corresponDIPg Jacobian matrix.
"""
function ∇ₚF(x, p, s)
    DIP, POP = unpackx(x)
    foo = zeros(promote_type(eltype(x), eltype(p)), n)
    if s == :umax
        foo[iPOP] .= ∂uptake_∂umax(DIP, p)
        foo[iDIP] .= -foo[iPOP]
    elseif s == :ku
        foo[iPOP] .= ∂uptake_∂ku(DIP, p)
        foo[iDIP] .= -foo[iPOP]
    elseif s == :w₀
        foo[iPOP] .= -S1 * POP
    elseif s == :w′
        foo[iPOP] .= -Sz * POP
    elseif s == :κ
        foo[iDIP] .= ∂remineralization_∂κ(POP, p)
        foo[iPOP] .= -foo[iDIP]
    elseif s == :DIPgeo
        foo[iDIP] .= ∂geores_∂DIPgeo(DIP, p)
    else
        error("There is no $s parameter")
    end
    return foo
end

