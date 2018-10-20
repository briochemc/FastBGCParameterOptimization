
@testset "Check cost functions output type" begin
    @testset "Check q!" begin
        @testset "with preprint" begin
            @testset "as a function of λ" begin
                @test q!(1λ₀)            isa Real
                @test q!(2λ₀ .+ 1e-50im) isa Complex
                @test q!(3λ₀ .+ ε)       isa Dual
            end
            @testset "as a function of p" begin
                @test q!(p₀) isa Real
            end
        end
        @testset "without preprint" begin
            @testset "as a function of λ" begin
                @test q!(4λ₀           , preprint = "  ") isa Real
                @test q!(5λ₀ .+ 1e-50im, preprint = "  ") isa Complex
                @test q!(6λ₀ .+ ε      , preprint = "  ") isa Dual
            end
            @testset "as a function of p" begin
                @test q!(p₀, preprint = "  ") isa Real
            end
        end
    end
    @testset "Check Dq!" begin
        @testset "with preprint" begin
            @test Dq!(1λ₀)            isa Array
            @test Dq!(2λ₀ .+ 1e-50im) isa Array
            @test Dq!(3λ₀ .+ ε)       isa Array
        end
        @testset "without preprint" begin
            @test Dq!(4λ₀           , preprint = "  ") isa Array
            @test Dq!(5λ₀ .+ 1e-50im, preprint = "  ") isa Array
            @test Dq!(6λ₀ .+ ε      , preprint = "  ") isa Array
        end
    end
    @testset "Check D2q!" begin
        @testset "with preprint" begin
            @test D2q!(1λ₀) isa Array
        end
        @testset "without preprint" begin
            @test D2q!(2λ₀, preprint = "  ") isa Array
        end
    end
end


