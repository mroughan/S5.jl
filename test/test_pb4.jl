@testset "IntermittentMapSymbols (PB4)" begin

    @testset "Constructor -- valid" begin
        g = IntermittentMapSymbols(1.6, [:a, :b], [0.3, 0.7]; burnin = 50)
        @test g.z == 1.6
        @test g.marginal == [0.3, 0.7]
        @test g.burnin == 50
    end

    @testset "Constructor -- argument errors" begin
        @test_throws ArgumentError IntermittentMapSymbols(1.0, [:a, :b])
        @test_throws ArgumentError IntermittentMapSymbols(1.6, [:a, :a])
        @test_throws ArgumentError IntermittentMapSymbols(1.6, [:a, :b], [0.4, 0.4])
        @test_throws ArgumentError IntermittentMapSymbols(1.6, [:a, :b], [-0.1, 1.1])
        @test_throws ArgumentError IntermittentMapSymbols(1.6, [:a, :b], [0.5, Inf])
        @test_throws ArgumentError IntermittentMapSymbols(1.6, [:a, :b]; burnin = -1)
    end

    @testset "generate -- output type, length, and finite-sample marginal" begin
        marginal = [0.2, 0.3, 0.5]
        g = IntermittentMapSymbols(1.6, [:a, :b, :c], marginal; burnin = 10)
        seq = generate(g, 1_003; rng = StableRNG(710))
        @test length(seq) == 1_003
        @test eltype(seq) == Symbol
        @test all(s in (:a, :b, :c) for s in seq)
        @test [count(==(s), seq) for s in (:a, :b, :c)] == bin_counts(marginal, 1_003)
    end

    @testset "generate -- rejects too-short sequences" begin
        g = IntermittentMapSymbols(1.6, [:a])
        @test_throws ArgumentError generate(g, 3)
    end

    @testset "target_marginal and capabilities" begin
        g = IntermittentMapSymbols(1.6, [:a, :b], [0.25, 0.75])
        @test target_marginal(g) == [0.25, 0.75]
        caps = control_capabilities(g)
        @test caps.alphabet == :exact
        @test caps.marginal == :finite_sample
        @test caps.lrd == :latent_empirical
    end

end
