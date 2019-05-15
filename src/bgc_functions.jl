
#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
T_DOP(p) = T_Circulation
const S₀ = buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_DOP, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(x, p)
    τg, xgeo = p.τg, p.xgeo
    return (xgeo .- x) / τg
end
# Uptake of phosphate (DIP)
soft_relu(x, p) = 0.5 * (tanh.(x / p.α) .+ 1) .* x
dsoft_relu(x, p) = -0.5 * (tanh.(x / p.α).^2 .- 1) .* x / p.α + 0.5tanh(x / p.α) + 0.5
function uptake(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺ = soft_relu(DIP, p)
    return Umax * DIP⁺ ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
# Remineralization DOP into DIP
function remineralization(DOP, p)
    κDOP = p.κDOP
    return κDOP * DOP
end
# Dissolution of POP into DOP
function dissolution(POP, p)
    κPOP = p.κPOP
    return κPOP * POP
end
# Add them up into sms functions (Sources Minus Sinks)
function sms_DIP(DIP, DOP, POP, p)
    return -uptake(DIP, p) + remineralization(DOP, p) + geores(DIP, p)
end
function sms_DOP(DIP, DOP, POP, p)
    σ = p.σ
    return σ * uptake(DIP, p) - remineralization(DOP, p) + dissolution(POP, p)
end
function sms_POP(DIP, DOP, POP, p)
    σ = p.σ
    return (1 - σ) * uptake(DIP, p) - dissolution(POP, p)
end
sms_all = (sms_DIP, sms_DOP, sms_POP) # bundles all the source-sink functions in a tuple
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
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺, dDIP⁺ = soft_relu(DIP), dsoft_relu(DIP)
    return sparse(Diagonal(Umax * ku * dDIP⁺ ./ (DIP⁺ .+ ku).^2 .* (z .≤ z₀)))
end
function ∂uptake_∂Umax(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺ = soft_relu(DIP)
    return DIP⁺ ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
function ∂uptake_∂ku(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺ = soft_relu(DIP)
    return -Umax * DIP⁺ ./ (DIP⁺ .+ ku).^2 .* (z .≤ z₀)
end

# Remineralization
function remineralizationJac(DOP, p)
    κDOP = p.κDOP
    return κDOP * I
end
function ∂remineralization_∂κDOP(DOP, p)
    return DOP
end
# Dissolution
function dissolutionJac(POP, p)
    κPOP = p.κPOP
    return κPOP * I
end
function ∂dissolution_∂κPOP(POP, p)
    return POP
end

# Indices for DIP, DOP, and POP (needed for Rate of change F(x,p))
const iDIP = 1:nb
const iDOP = iDIP .+ nb
const iPOP = iDOP .+ nb
function unpackx(x)
    DIP = x[iDIP]
    DOP = x[iDOP]
    POP = x[iPOP]
    return DIP, DOP, POP
end
# add method to deal with arrays for numJac
#function unpackx(x::Array{<:Number,2})
#    DIP = x[iDIP,:]
#    POP = x[iPOP,:]
#    return DIP, POP
#end

# F and ∇ₓF by hand
function F(x, p)
    DIP, DOP, POP = unpackx(x)
    σ = p.σ
    u = uptake(DIP, p)
    r = remineralization(DOP, p)
    d = dissolution(POP, p)
    return [ -T_DIP(p) * DIP         - u     + r + geores(DIP, p) ;
             -T_DOP(p) * DOP +    σ  * u + d - r                  ;
             -T_POP(p) * POP + (1-σ) * u - d                      ]
end

function ∇ₓF(x, p)
    DIP, DOP, POP = unpackx(x)
    σ = p.σ
    uJac = uptakeJac(DIP, p)
    rJac = remineralizationJac(DOP, p)
    dJac = dissolutionJac(POP, p)
    foo = [ -T_DIP(p) - uJac + georesJac(p)              rJac                   ;
                   σ  * uJac                 -T_DOP(p) - rJac              dJac ;
              (1 - σ) * uJac                                   -T_POP(p) - dJac ]
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
    DIP, DOP, POP = unpackx(x)
    σ = p.σ
    ∇ₚF = zeros(promote_type(eltype(x), eltype(p)), n)
    if s == :Umax
        ∇ₚF[iDIP] .= -∂uptake_∂Umax(DIP, p)
        ∇ₚF[iDOP] .=   -σ   * ∇ₚF[iDIP]
        ∇ₚF[iPOP] .= -(1-σ) * ∇ₚF[iDIP]
    elseif s == :ku
        ∇ₚF[iDIP] .= -∂uptake_∂ku(DIP, p)
        ∇ₚF[iDOP] .=   -σ   * ∇ₚF[iDIP]
        ∇ₚF[iPOP] .= -(1-σ) * ∇ₚF[iDIP]
    elseif s == :w₀
        ∇ₚF[iPOP] .= -S₀ * POP
    elseif s == :w′
        ∇ₚF[iPOP] .= -S′ * POP
    elseif s == :κDOP
        ∇ₚF[iDIP] .= ∂remineralization_∂κDOP(DOP, p)
        ∇ₚF[iDOP] .= -∇ₚF[iDIP]
    elseif s == :κPOP
        ∇ₚF[iDOP] .= ∂dissolution_∂κPOP(POP, p)
        ∇ₚF[iPOP] .= -∇ₚF[iDOP]
    elseif s == :xgeo
        ∇ₚF[iDIP] .= ∂geores_∂xgeo(DIP, p)
    else
        error("There is no $s parameter")
    end
    return ∇ₚF
end

