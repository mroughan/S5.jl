@testset "FSS (MB3)" begin

    @testset "Constructor — valid inputs" begin
        g = FSS(1.5, [:a, :b, :c])
        @test g.alpha == 1.5
        @test isapprox((3 - g.alpha) / 2, 0.75)   # H = 0.75
        @test g.x_min == 1.0
        @test length(g.rates) == 3
        @test all(==(1.0), g.rates)

        g2 = FSS(1.2, [:x, :y]; rates = [1.0, 3.0], x_min = 0.5)
        @test g2.rates == [1.0, 3.0]
        @test g2.x_min == 0.5
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError FSS(1.0, [:a, :b])                        # alpha = 1
        @test_throws ArgumentError FSS(2.0, [:a, :b])                        # alpha = 2
        @test_throws ArgumentError FSS(1.5, [:a, :b]; rates = [-1.0, 1.0])  # negative rate
        @test_throws ArgumentError FSS(1.5, [:a, :b, :c]; rates = [1.0, 1.0]) # length mismatch
        @test_throws ArgumentError FSS(1.5, [:a]; x_min = 0.0)              # x_min = 0
    end

    @testset "generate — output type and length" begin
        g   = FSS(1.5, ['x', 'y', 'z'])
        seq = generate(g, 2_000; rng = rng)
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('x', 'y', 'z') for c in seq)
    end

    @testset "generate — uniform marginal" begin
        g   = FSS(1.5, [:a, :b, :c])
        seq = generate(g, 9_000; rng = rng)
        for s in (:a, :b, :c)
            @test isapprox(count(==(s), seq) / 9_000, 1/3; atol = 0.04)
        end
    end

    @testset "generate — rate-controlled marginal" begin
        # rates = [1, 3] → symbol :b should appear ~75 % of the time
        g   = FSS(1.5, [:a, :b]; rates = [1.0, 3.0])
        seq = generate(g, 20_000; rng = rng)
        @test isapprox(count(==(:a), seq) / 20_000, 0.25; atol = 0.03)
        @test isapprox(count(==(:b), seq) / 20_000, 0.75; atol = 0.03)
    end

    @testset "generate — rejects n < 1" begin
        g = FSS(1.5, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "Pareto sampler — tail behaviour (statistical)" begin
        # _pareto_sample should produce heavy-tailed draws; check empirical mean
        # and that a fraction of samples exceeds the theoretical 90th percentile.
        alpha = 1.5
        x_min = 1.0
        N     = 50_000
        draws = [S5._pareto_sample(MersenneTwister(i), alpha, x_min) for i in 1:N]
        @test all(>(x_min), draws)                  # all draws > x_min
        p90   = x_min / 0.1^(1 / alpha)            # theoretical 90th percentile
        @test isapprox(count(>(p90), draws) / N, 0.1; atol = 0.02)
    end

end
