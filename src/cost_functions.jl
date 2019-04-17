
#===================================
Generate 𝑓 and ∇ₓ𝑓
===================================#
# hyper parameters
ωPO4, ωPOP = 1.0, 0.0   # no cost for POP
ωs = (ωPO4, ωPOP)       # tracers weight
ωp = 1e-4               # parameter weight
# mean and variance of PO4 and POP obs
xobs = (PO4obs, PO4obs)         # observations for tracers
σ²xobs = (σ²PO4obs, σ²PO4obs)   # variance of observations for tracers
# Build 𝑓 and ∇ₓ𝑓
f, ∇ₓf = TransportMatrixTools.mismatch_function_and_Jacobian(ωs, xobs, σ²xobs, v, ωp, logpobs, σ²logpobs)

