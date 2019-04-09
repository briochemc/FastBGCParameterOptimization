#=============================================
    Test objective, gradient, and Hessian
    for each method
=============================================#

list_∇ᵏf̂ = [:f̂, :∇f̂, :∇²f̂]
list_methods = [:F1, :DUAL, :CSD, :FD1, :HYPER, :FD2]
λ = exp(randn()) .* λ₀

@testset "Testing F1 method against" begin
    @testset "$m method" for m in list_methods[2:end] # Test DUAL, CSD, FD1
        @testset "$∇ᵏf̂" for ∇ᵏf̂ in list_∇ᵏf̂
            m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
            F1_∇ᵏf̂ = Symbol("F1_" * string(∇ᵏf̂))
            rtol = m ∈ list_methods[[4,6]] ? 1e-4 : 1e-9
            @eval @test isapprox($m_∇ᵏf̂(λ), $F1_∇ᵏf̂(λ), rtol=$rtol)
        end
    end
end

#=============================================
    Table of relative mismatches
    between each method and the F1 method
=============================================#

@testset "Print relative mismatch of F1 method vs" begin
    str = ""
    relative_difference(x, y) = norm(x - y) / norm(x)
    for ∇ᵏf̂ in list_∇ᵏf̂
        str = str * "    " * string(∇ᵏf̂) * "\n"
        for m in list_methods[2:end] # Test DUAL, CSD, FD1
            m_∇ᵏf̂ = Symbol(string(m) * "_" * string(∇ᵏf̂))
            F1_∇ᵏf̂ = Symbol("F1_" * string(∇ᵏf̂))
            str = str * "        " * (@sprintf "%5s" string(m)) * " method: "
            @eval foo = relative_difference($m_∇ᵏf̂(λ), $F1_∇ᵏf̂(λ))
            str = str * (@sprintf "%8.2e" foo) * "\n"
        end
    end
    println(str)
end



