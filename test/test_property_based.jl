@testset "Composable property-based generators" begin

    @testset "source and symbolizer metadata" begin
        q = QuantileSymbolizer([:a, :b], [0.25, 0.75])
        a = ArgmaxSymbolizer([:a, :b, :c])
        P1 = [0.9 0.1; 0.2 0.8]
        P2 = [0.3 0.7; 0.6 0.4]
        m = MarkovRegimeSymbolizer([:a, :b], [P1, P2])

        @test q isa Symbolizer
        @test a isa Symbolizer
        @test m isa Symbolizer
        @test SpectralFGNSource(0.8) isa LatentSource
        @test HaarLRDSource(0.8) isa LatentSource
        @test IntermittentMapSource(1.6) isa LatentSource

        @test latent_width(q) == 1
        @test latent_width(a) == 3
        @test latent_width(m) == 1
        @test target_marginal(PropertyBasedGenerator(SpectralFGNSource(0.8), q)) ==
              [0.25, 0.75]
    end

    @testset "constructor errors" begin
        @test_throws ArgumentError SpectralFGNSource(0.5)
        @test_throws ArgumentError HaarLRDSource(1.0)
        @test_throws ArgumentError IntermittentMapSource(1.0)
        @test_throws ArgumentError QuantileSymbolizer([:a, :a])
        @test_throws ArgumentError ArgmaxSymbolizer([:a, :b], [0.4, 0.4])
        @test_throws ArgumentError MarkovRegimeSymbolizer(
            [:a, :b], [[0.5 0.5; 0.5 0.5]]; regime_weights = [0.4, 0.4])
        @test_throws ArgumentError PropertyBasedGenerator(
            IntermittentMapSource(1.6), ArgmaxSymbolizer([:a, :b]))
    end

    @testset "spectral source with quantile symbolizer" begin
        marginal = [0.2, 0.3, 0.5]
        g = PropertyBasedGenerator(SpectralFGNSource(0.8),
                                   QuantileSymbolizer([:a, :b, :c], marginal))
        seq = generate(g, 1003; rng = StableRNG(910))
        @test length(seq) == 1003
        @test eltype(seq) == Symbol
        @test [count(==(s), seq) for s in (:a, :b, :c)] ==
              bin_counts(marginal, length(seq))
        @test generate(g, 128; rng = StableRNG(911)) ==
              generate(g, 128; rng = StableRNG(911))
    end

    @testset "spectral source with argmax symbolizer" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        g = PropertyBasedGenerator(
            SpectralFGNSource(0.8),
            ArgmaxSymbolizer(alphabet, marginal; calibration_iters = 35))
        seq = generate(g, 10_000; rng = StableRNG(912))
        @test length(seq) == 10_000
        @test all(s in alphabet for s in seq)
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.04
    end

    @testset "latent generation and direct symbolization" begin
        latent = generate_latent(SpectralFGNSource(0.75), 16, 2;
                                 rng = StableRNG(913))
        @test size(latent) == (2, 16)

        symbols = symbolize(QuantileSymbolizer([:lo, :hi]),
                            reshape(collect(1.0:8.0), 1, 8);
                            rng = StableRNG(914))
        @test symbols == [:lo, :lo, :lo, :lo, :hi, :hi, :hi, :hi]
    end

    @testset "Markov regime symbolizer" begin
        P1 = [0.9 0.1; 0.2 0.8]
        P2 = [0.3 0.7; 0.6 0.4]
        g = PropertyBasedGenerator(
            HaarLRDSource(0.8; cascade_depth = 3),
            MarkovRegimeSymbolizer([:a, :b], [P1, P2]; regime_weights = [0.4, 0.6]))
        seq = generate(g, 256; rng = StableRNG(915))
        @test length(seq) == 256
        @test all(s in (:a, :b) for s in seq)
        @test target_marginal(g) ≈ 0.4 .* stationary_distribution(P1) .+
                                  0.6 .* stationary_distribution(P2)
        @test control_capabilities(g).bigram == :per_regime
    end

    @testset "intermittent source with quantile symbolizer" begin
        marginal = [0.4, 0.6]
        g = PropertyBasedGenerator(
            IntermittentMapSource(1.6; burnin = 10),
            QuantileSymbolizer(['A', 'B'], marginal))
        seq = generate(g, 1001; rng = StableRNG(916))
        @test length(seq) == 1001
        @test eltype(seq) == Char
        @test [count(==(s), seq) for s in ('A', 'B')] ==
              bin_counts(marginal, length(seq))
    end

end
