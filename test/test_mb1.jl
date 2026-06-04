# Simple lag-k ACF for an indicator series; no external dependencies.
function _acf(x::AbstractVector{<:Real}, k::Int)
    n = length(x)
    μ = sum(x) / n
    v = sum((xi - μ)^2 for xi in x) / n
    v == 0.0 && return 0.0
    return sum((x[i] - μ) * (x[i + k] - μ) for i in 1:(n - k)) / (n * v)
end

@testset "LAMP (MB1)" begin

    @testset "Constructor — valid inputs" begin
        g = LAMP(0.5, [:a, :b, :c])
        @test g.beta == 0.5
        @test g.d    == 1000
        @test isapprox(sum(g.marginal), 1.0)
        @test isapprox(sum(g.weights), 1.0)
        @test g.weights[1] > g.weights[2] > g.weights[end]   # decreasing

        g2 = LAMP(0.3, [:x, :y]; d = 200)
        @test g2.d == 200
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError LAMP(0.0,  [:a, :b])              # beta = 0
        @test_throws ArgumentError LAMP(1.0,  [:a, :b])              # beta = 1
        @test_throws ArgumentError LAMP(0.5,  [:a, :b], [0.3, 0.5]) # sum ≠ 1
        @test_throws ArgumentError LAMP(0.5,  [:a, :b, :c], [0.5, 0.5]) # length mismatch
        @test_throws ArgumentError LAMP(0.5,  [:a, :b]; d = 0)      # d < 1
    end

    @testset "generate — output type and length" begin
        g   = LAMP(0.5, ['a', 'b', 'c']; d = 100)
        seq = generate(g, 2_000; rng = rng)
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b', 'c') for c in seq)
    end

    @testset "generate — marginal frequencies (uniform)" begin
        g   = LAMP(0.5, [:a, :b, :c]; d = 200)
        seq = generate(g, 9_000; rng = rng)
        for s in (:a, :b, :c)
            @test isapprox(count(==(s), seq) / 9_000, 1/3; atol = 0.04)
        end
    end

    @testset "generate — marginal frequencies (non-uniform)" begin
        g   = LAMP(0.4, [1, 2], [0.3, 0.7]; d = 100)
        seq = generate(g, 6_000; rng = rng)
        @test isapprox(count(==(1), seq) / 6_000, 0.3; atol = 0.04)
        @test isapprox(count(==(2), seq) / 6_000, 0.7; atol = 0.04)
    end

    @testset "generate — rejects n < 1" begin
        g = LAMP(0.5, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "ACF is positive and decreasing (statistical)" begin
        # A sequence with strong LRD (beta=0.3) should show positive, decreasing ACF.
        g   = LAMP(0.3, [:a, :b]; d = 500)
        seq = generate(g, 15_000; rng = MersenneTwister(200))
        x   = Float64.(seq .== :a)

        acf1  = _acf(x, 1)
        acf10 = _acf(x, 10)
        acf50 = _acf(x, 50)

        @test acf1  > 0
        @test acf10 > 0
        @test acf1  > acf10 > acf50
    end

end
