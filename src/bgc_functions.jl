
#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
const S₀ = buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(x, p)
    τg, xgeo = p.τg, p.xgeo
    return (xgeo .- x) / τg
end
# Uptake of phosphate (DIP)
relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return 1/τu * DIP⁺.^2 ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
# Remineralization DOP into DIP
function remineralization(POP, p)
    κ = p.κ
    return κ * POP
end
# Add them up into sms functions (Sources Minus Sinks)
sms_DIP(DIP, POP, p) = -uptake(DIP, p) + remineralization(POP, p) + geores(DIP, p)
sms_POP(DIP, POP, p) =  uptake(DIP, p) - remineralization(POP, p)
sms_all = (sms_DIP, sms_POP) # bundles all the source-sink functions in a tuple
 # generates the state function (and its Jacobian!)

#===========================================
AIBECS F and ∇ₓF
===========================================#
A_F, A_∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

#===========================================
Analytical derivatives
Required for CSD, DUAL, and HYPER methods
===========================================#

# Geological restoring
function georesJac(p)
    τg = p.τg
    return -I / τg
end
function ∂geores_∂xgeo(x, p)
    τg = p.τg
    return 1 / τg
end

# Uptake
drelu(x) = (x .≥ 0)
function uptakeJac(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺, dDIP⁺ = relu(DIP), drelu(DIP)
    return sparse(Diagonal(1/τu * dDIP⁺ .* DIP⁺ .* (DIP⁺ .+ 2ku) ./ (DIP⁺ .+ ku).^2 .* (z .≤ z₀)))
end
function ∂uptake_∂τu(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return -1/τu^2 * DIP⁺.^2 ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
function ∂uptake_∂ku(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return -1/τu * DIP⁺.^2 ./ (DIP⁺ .+ ku).^2 .* (z .≤ z₀)
end

# Remineralization
function remineralizationJac(POP, p)
    κ = p.κ
    return κ * I
end
function ∂remineralization_∂κ(POP, p)
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
#function unpackx(x::Array{<:Number,2})
#    DIP = x[iDIP,:]
#    POP = x[iPOP,:]
#    return DIP, POP
#end

# F and ∇ₓF by hand
function F(x, p)
    DIP, POP = unpackx(x)
    u = uptake(DIP, p)
    r = remineralization(POP, p)
    return [ -T_DIP(p) * DIP - u + r + geores(DIP, p) ;
             -T_POP(p) * POP + u - r                  ]
end

function ∇ₓF(x, p)
    DIP, POP = unpackx(x)
    uJac = uptakeJac(DIP, p)
    rJac = remineralizationJac(POP, p)
    foo = [ -T_DIP(p) - uJac + georesJac(p)              rJac  ;
                        uJac                 -T_POP(p) - rJac ]
    return foo
end

using Flatten
"""
    ∇ₚF(x, p)

Evaluates the jacobian of `f` with respect to `p`.
Concatenates `∇ₚF(x, p::Para, s::Symbol)` for all optimizable parameter symbols `s`.
"""
∇ₚF(x, p) = reduce(hcat, ∇ₚF(x, p, s) for s in fieldnameflatten(p))

"""
    ∇ₚF(x, p, s::Symbol)

Evaluates the derivative of `F` with respect to `p.s` where `s` is the name of the parameter (of type `Symbol`).

You should fill this function with all the first derivatives of f with respoect to each parameter.
Called without the symbol `s`, this function will loop through all the optimizable parameters and create the corresponDIPg Jacobian matrix.
"""
function ∇ₚF(x, p, s)
    DIP, POP = unpackx(x)
    ∇ₚF = zeros(promote_type(eltype(x), eltype(p)), n)
    if s == :τu
        ∇ₚF[iDIP] .= -∂uptake_∂τu(DIP, p)
        ∇ₚF[iPOP] .= -∇ₚF[iDIP]
    elseif s == :ku
        ∇ₚF[iDIP] .= -∂uptake_∂ku(DIP, p)
        ∇ₚF[iPOP] .= -∇ₚF[iDIP]
    elseif s == :w₀
        ∇ₚF[iPOP] .= -S₀ * POP
    elseif s == :w′
        ∇ₚF[iPOP] .= -S′ * POP
    elseif s == :κ
        ∇ₚF[iDIP] .= ∂remineralization_∂κ(POP, p)
        ∇ₚF[iPOP] .= -∇ₚF[iDIP]
    elseif s == :xgeo
        ∇ₚF[iDIP] .= ∂geores_∂xgeo(DIP, p)
    else
        error("There is no $s parameter")
    end
    return ∇ₚF
end

