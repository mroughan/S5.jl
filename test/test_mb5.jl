@testset "CalibratedAdditiveMarkov (MB1c)" begin

    @testset "Constructor -- valid" begin
        g = CalibratedAdditiveMarkov(0.4, [:a, :b], [0.3, 0.7];
                                     d = 25, strength = 0.6)
        @test g.beta == 0.4
        @test g.d == 25
        @test g.strength == 0.6
        @test isapprox(sum(g.weights), 1.0)
    end

    @testset "Constructor -- argument errors" begin
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.0, [:a, :b])
        @test_throws ArgumentError CalibratedAdditiveMarkov(1.0, [:a, :b])
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :a])
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :b], [0.4, 0.4])
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :b], [0.5, NaN])
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :b]; d = 0)
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :b]; strength = -0.1)
        @test_throws ArgumentError CalibratedAdditiveMarkov(0.5, [:a, :b]; strength = 1.1)
    end

    @testset "generate -- output and iid strength zero" begin
        g = CalibratedAdditiveMarkov(0.5, ['x', 'y'], [0.2, 0.8];
                                     d = 20, strength = 0.0)
        seq = generate(g, 10_000; rng = StableRNG(720))
        @test length(seq) == 10_000
        @test eltype(seq) == Char
        @test all(c in ('x', 'y') for c in seq)
        @test total_variation(empirical_marginal(seq, ['x', 'y']), [0.2, 0.8]) < 0.02
    end

    @testset "history strength increases repeats" begin
        iid = CalibratedAdditiveMarkov(0.5, [:a, :b]; d = 40, strength = 0.0)
        persistent = CalibratedAdditiveMarkov(0.5, [:a, :b]; d = 40,
                                              strength = 0.9)
        seq_iid = generate(iid, 5_000; rng = StableRNG(721))
        seq_persistent = generate(persistent, 5_000; rng = StableRNG(721))
        repeat_rate(seq) = count(seq[i] == seq[i - 1] for i in 2:length(seq)) /
                           (length(seq) - 1)
        @test repeat_rate(seq_persistent) > repeat_rate(seq_iid) + 0.05
    end

    @testset "target_marginal and capabilities" begin
        g = CalibratedAdditiveMarkov(0.5, [:a, :b], [0.25, 0.75])
        @test target_marginal(g) == [0.25, 0.75]
        caps = control_capabilities(g)
        @test caps.alphabet == :exact
        @test caps.marginal == :centered_target
        @test caps.lrd == :finite_history
    end

end

@testset "DuplicationMutation (MB5)" begin

    function indicator_acf(seq, symbol, lag)
        x = Float64.(seq .== symbol)
        μ = sum(x) / length(x)
        v = sum((xi - μ)^2 for xi in x) / length(x)
        return sum((x[i] - μ) * (x[i + lag] - μ)
                   for i in 1:(length(x) - lag)) / ((length(x) - lag) * v)
    end

    @testset "Constructor -- valid" begin
        g = DuplicationMutation(1.5, ['A', 'C', 'G', 'T'];
                                mutation_probability = 0.02,
                                seed_length = 12,
                                max_block_length = 50)
        @test g.alpha == 1.5
        @test g.mutation_probability == 0.02
        @test g.seed_length == 12
        @test g.max_block_length == 50
    end

    @testset "Constructor -- argument errors" begin
        @test_throws ArgumentError DuplicationMutation(1.0, [:a, :b])
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :a])
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b], [0.4, 0.4])
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b], [0.5, Inf])
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b];
                                                       mutation_probability = -0.1)
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b];
                                                       mutation_probability = 1.1)
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b]; seed_length = 0)
        @test_throws ArgumentError DuplicationMutation(1.5, [:a, :b]; max_block_length = 0)
    end

    @testset "generate -- output type and length" begin
        g = DuplicationMutation(1.5, ['A', 'C', 'G', 'T'];
                                mutation_probability = 0.05,
                                seed_length = 16,
                                max_block_length = 100)
        seq = generate(g, 2_000; rng = StableRNG(730))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c in ('A', 'C', 'G', 'T') for c in seq)
    end

    @testset "full mutation recovers marginal approximately" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        g = DuplicationMutation(1.5, alphabet, marginal;
                                mutation_probability = 1.0,
                                seed_length = 8,
                                max_block_length = 50)
        seq = generate(g, 20_000; rng = StableRNG(731))
        @test total_variation(empirical_marginal(seq, alphabet), marginal) < 0.02
    end

    @testset "copy distance creates decaying dependence" begin
        g = DuplicationMutation(1.4, [:a, :b];
                                mutation_probability = 0.02,
                                seed_length = 64,
                                max_block_length = 5_000)
        seq = generate(g, 30_000; rng = StableRNG(732))
        near = indicator_acf(seq, :a, 10)
        middle = indicator_acf(seq, :a, 200)
        far = indicator_acf(seq, :a, 2_000)
        @test near > middle
        @test middle > far
        @test near > 0.05
    end

    @testset "target_marginal and capabilities" begin
        g = DuplicationMutation(1.5, [:a, :b], [0.25, 0.75])
        @test target_marginal(g) == [0.25, 0.75]
        caps = control_capabilities(g)
        @test caps.alphabet == :exact
        @test caps.marginal == :mutation_target
        @test caps.lrd == :empirical
    end

end
