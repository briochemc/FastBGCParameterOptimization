#===================================
Generic model has the form
    ∂𝒙/∂𝑡 = 𝑭(𝒙,𝒑) = 0,
where
    𝑭(𝒙,𝒑) = 𝐓(𝒑) 𝒙 + 𝑮(𝒙,𝒑).
- 𝐓(𝒑) is a block-diagonal matrix of the transport matrices for each tracer.
    It is constructed from a tuple of functions of 𝒑 that the user must supply.
- 𝑮(𝒙,𝒑) is the nonlinear local sources minus sinks for each tracer.
    It is constructed from a tuple of functions of (𝒙,𝒑) that the user must supply.
===================================#

# Load useful OCIM constants
const nb, DIV, Iabove, v, z, ztop = constants(wet3d, grd)

#===================================
Geological restoring for PO4
===================================#
# Load WOA mean and variance for PO4
const PO4obs, σ²PO4obs = TransportMatrixTools.build_μ_and_σ²_from_WOA(wet3d, grd, "PO4")
# Compute mean observed PO4 (averaged over entire ocean volume)
const meanPO4obs = TransportMatrixTools.vmean(PO4obs, v, findall(isfinite.(σ²PO4obs)))
# Geological restoring
function geores(x, p)
    τg = p.τg
    return (meanPO4obs .- x) / τg
end

#===================================
Uptake of PO4
===================================#
# postivie part
relu(x) = (x .≥ 0) .* x
# Michaelis-Menten function
mm(x, μ, k) = μ * x ./ (x .+ k)
# Depth of the base of the euphotic zone
const z₀ = 85 # 𝑧₀ = 85m ⟹  2 layers
# Uptake
function uptake(x, p)
    umax, ku = p.umax, p.ku
    return (z .≥ z₀) .* mm(relu(x), umax, ku)
end

#===================================
Remineralization of POP
===================================#
function remineralization(x, p)
    κ = p.κ
    return κ * x
end

#===================================
Transport matrix
for PO4
===================================#
T_PO4(p) = T_OCIM

#===================================
Transport matrix
for POP
===================================#
# PFD transport (needed for Rate of change F(x,p))
const S₀= buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end

#===================================
Components of 𝑮(𝒙,𝒑), which is
the nonlinear part of 𝑭(𝒙,𝒑)
===================================#
G_PO4(PO4, POP, p) = -uptake(PO4, p) + remineralization(POP, p) + geores(PO4, p)
G_POP(PO4, POP, p) =  uptake(PO4, p) - remineralization(POP, p)

Ts = (T_PO4, T_POP)
Gs = (G_PO4, G_POP)

# number of tracers
const nt = length(Ts)

#===================================
Generate 𝑭 and ∇ₓ𝑭
===================================#
F, ∇ₓF = TransportMatrixTools.multiTracer.build_F_and_∇ₓF(Ts, Gs, nt, nb)

#===================================
Generate 𝑓 and ∇ₓ𝑓
===================================#
# hyper parameters
ωPO4, ωPOP = 1.0, 0.0   # no cost for POP
ωs = (ωPO4, ωPOP)       # tracers weight
ωp = 1e-4               # parameter weight
xobs = (PO4obs, PO4obs)         # observations for tracers
σ²xobs = (σ²PO4obs, σ²PO4obs)   # variance of observations for tracers
# TODO pobs = ??
# TODO σ²p = ??
# TODO f, ∇ₓf = build_f_and_∇ₓf(ωs, xobs, σ²xobs, v, ωp, pobs, σ²pobs)
