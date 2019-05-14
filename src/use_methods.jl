

#=============================================
    Assign function name to method
=============================================#

objective_gradient_hessian = [:objective, :gradient, :hessian]
list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:AF1, :F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
Para = AIBECS.Parameters

# Objective
AF1_f̂(p::Para; preprint=" ") = F1.objective(f, F, A_∇ₓF, AF1_mem, p, CTKAlg(), nrm=nrm, preprint=preprint)
F1_f̂(p::Para; preprint=" ") = F1.objective(f, F, ∇ₓF, F1_mem, p, CTKAlg(), nrm=nrm, preprint=preprint)
for m in list_methods[3:end] # All methods except F1
    m_f̂ = Symbol(string(m) * "_f̂")
    msol = Symbol(string(m) * "Sol")
    @eval $m_f̂(p::Para; preprint=" ") = $m.objective(f, F, ∇ₓF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
end

# Gradient and Hessian
for (∇ᵏf̂, f_g_h) in zip(list_∇ᵏf̂[2:3], objective_gradient_hessian[2:3])
    # AIBECS + F1 method
    m_∇ᵏf̂ = Symbol("AF1_" * string(∇ᵏf̂))
    @eval $m_∇ᵏf̂(p::Para; preprint=" ") = F1.$f_g_h(f, F, ∇ₓf, A_∇ₓF, AF1_mem, p, CTKAlg(), nrm=nrm, preprint=preprint)
    # F1 method
    m_∇ᵏf̂ = Symbol("F1_" * string(∇ᵏf̂))
    @eval $m_∇ᵏf̂(p::Para; preprint=" ") = F1.$f_g_h(f, F, ∇ₓf, ∇ₓF, F1_mem, p, CTKAlg(), nrm=nrm, preprint=preprint)
    # DUAL, CSD, and FD1 methods
    for m in list_methods[3:5]
        msol = Symbol(string(m) * "Sol")
        m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
        @eval $m_∇ᵏf̂(p::Para; preprint=" ") = $m.$f_g_h(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
    end
    # HYPER and FD2 methods
    for m in list_methods[6:7]
        msol = Symbol(string(m) * "Sol")
        m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
        @eval $m_∇ᵏf̂(p::Para; preprint=" ") = $m.$f_g_h(f, F, ∇ₓF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
    end
end

#=============================================
    Change of variable from p to λ
=============================================#

# Change variable from vector λ to Parameters p
λ2p(λ) = AIBECS.opt_para(exp.(λ))
# But derivatives of this change of variables
# must outputs vectors!
∇λ2p(λ) = exp.(λ)
∇²λ2p(λ) = exp.(λ)
# Reverse change of variables from Parameters p to vector λ
p2λ(p) = log.(optvec(p))
const λ₀ = p2λ(p₀)
# Then apply the composition rule to get the functions of λ
for m in list_methods
    m_f̂ = Symbol(string(m) * "_f̂")
    m_∇f̂ = Symbol(string(m) * "_∇f̂")
    m_∇²f̂ = Symbol(string(m) * "_∇²f̂")
    @eval begin
        $m_f̂(λ; preprint=" ") = $m_f̂(λ2p(λ), preprint=preprint)
        $m_∇f̂(λ; preprint=" ") = $m_∇f̂(λ2p(λ), preprint=preprint) * Diagonal(∇λ2p(λ))
        function $m_∇²f̂(λ; preprint=" ")
            ∇p = Diagonal(∇λ2p(λ)) # for variable change
            ∇²p = Diagonal(∇²λ2p(λ)) # for variable change
            G = vec($m_∇f̂(λ2p(λ), preprint=preprint))
            H = $m_∇²f̂(λ2p(λ), preprint=preprint)
            return ∇p * H * ∇p + Diagonal(G) * ∇²p
        end
    end
end

#=============================================
    Objective and derivatives with storage
=============================================#

for m in list_methods
    m_∇f̂ = Symbol(string(m) * "_∇f̂")
    m_∇²f̂ = Symbol(string(m) * "_∇²f̂")
    @eval $m_∇f̂(s, λ) = s[1:m] .= vec($m_∇f̂(λ))
    @eval $m_∇²f̂(s, λ) = s[1:m, 1:m] .= $m_∇²f̂(λ)
end

#=============================================
    Extra requirements for methods to run
=============================================#

# Because CSD will convert p.z₀ to the Complex type
# Solution: Assuming ℜ is defined, I overload isless
Base.isless(a::Real, b::Complex{Float64}) = isless(a, ℜ(b))
