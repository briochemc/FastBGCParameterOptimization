

#=============================================
    Assign function name to method
=============================================#

list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]

# Objective
F1_f̂(p::Para; preprint=" ") = F1.f̂(f, F, ∇ₓF, F1buf, p, CTKAlg(), nrm=nrm, preprint=preprint)
for m in list_methods[2:end] # All methods except F1
    m_f̂ = Symbol(string(m) * "_f̂")
    msol = Symbol(string(m) * "Sol")
    @eval $m_f̂(p::Para; preprint=" ") = $m.f̂(f, F, ∇ₓF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
end

# Gradient and Hessian
for ∇ᵏf̂ in list_∇ᵏf̂[2:3]
    m_∇ᵏf̂ = Symbol("F1_" * string(∇ᵏf̂))
    @eval  $m_∇ᵏf̂(p::Para; preprint=" ") = F1.$∇ᵏf̂(f, F, ∇ₓf, ∇ₓF, F1buf, p, CTKAlg(), nrm=nrm, preprint=preprint)
    for m in list_methods[2:4] # DUAL, CSD, and FD1
        msol = Symbol(string(m) * "Sol")
        m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
        @eval $m_∇ᵏf̂(p::Para; preprint=" ") = $m.$∇ᵏf̂(f, F, ∇ₓf, ∇ₓF, ∇ₚf, ∇ₚF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
    end
    for m in list_methods[5:6] # HYPER and FD2
        msol = Symbol(string(m) * "Sol")
        m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
        @eval $m_∇ᵏf̂(p::Para; preprint=" ") = $m.$∇ᵏf̂(f, F, ∇ₓF, $msol, p, CTKAlg(), nrm=nrm, preprint=preprint)
    end
end

#=============================================
    Change of variable from p to λ
=============================================#

for m in list_methods
    m_f̂ = Symbol(string(m) * "_f̂")
    m_∇f̂ = Symbol(string(m) * "_∇f̂")
    m_∇²f̂ = Symbol(string(m) * "_∇²f̂")
    @eval begin
        $m_f̂(λ; preprint=" ") = $m_f̂(λ2p(λ), preprint=preprint)
        $m_∇f̂(λ; preprint=" ") = $m_∇f̂(λ2p(λ), preprint=preprint) * diagm(0 => ∇λ2p(λ))
        function $m_∇²f̂(λ; preprint=" ")
            ∇p = diagm(0 => ∇λ2p(λ)) # for variable change
            ∇²p = diagm(0 => ∇²λ2p(λ)) # for variable change
            H = $m_∇²f̂(λ2p(λ), preprint=preprint)
            G = vec($m_∇f̂(λ2p(λ), preprint=preprint))
            return ∇p * H * ∇p + diagm(0 => G) * ∇²p
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



