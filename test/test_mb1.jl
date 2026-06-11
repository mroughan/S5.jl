function _acf(x::AbstractVector{<:Real}, k::Int)
    n = length(x)
    μ = sum(x) / n
    v = sum((xi - μ)^2 for xi in x) / n
    v == 0.0 && return 0.0
    return sum((x[i] - μ) * (x[i + k] - μ) for i in 1:(n - k)) / (n * v)
end

@testset "DyadicLAMP (MB1b)" begin

    @testset "Constructor — valid" begin
        g = DyadicLAMP(0.5, [:a, :b, :c]; d = 10_000)
        @test g.beta == 0.5
        @test g.d == 10_000
        @test g.epsilon == 0.01
        @test isapprox(sum(g.weights), 1.0)
        @test g.prefix_weights[1] == 0.0
        @test g.prefix_weights[end] == 1.0
        @test g.weights[1] > g.weights[2] > g.weights[end]
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError DyadicLAMP(0.0, [:a, :b])
        @test_throws ArgumentError DyadicLAMP(1.0, [:a, :b])
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :b], [0.3, 0.5])
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :b]; d = 0)
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :b]; epsilon = -0.1)
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :b]; epsilon = 1.1)
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :a])
        @test_throws ArgumentError DyadicLAMP(0.5, [:a, :b];
                                             transition_matrix = [0.5 0.5])
    end

    @testset "generate — output type and length" begin
        g = DyadicLAMP(0.5, ['a', 'b', 'c']; d = 10_000)
        seq = generate(g, 2_000; rng = StableRNG(11))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b', 'c') for c in seq)
    end

    @testset "generate — supports d greater than n without random prehistory" begin
        g = DyadicLAMP(0.5, [:a, :b], [1.0, 0.0]; d = 1_000_000, epsilon = 0.0)
        seq = generate(g, 20; rng = StableRNG(124))
        @test seq == fill(:a, 20)
    end

    @testset "generate — transition matrix controls repeated symbols" begin
        n = 2_000
        sticky = DyadicLAMP(0.5, [:a, :b], [0.5, 0.5];
                            d = 10_000,
                            epsilon = 0.0,
                            transition_matrix = lamp_repeat_transition([0.5, 0.5];
                                                                       repeat_probability = 0.95))
        iidlike = DyadicLAMP(0.5, [:a, :b], [0.5, 0.5];
                             d = 10_000,
                             epsilon = 0.0,
                             transition_matrix = lamp_repeat_transition([0.5, 0.5];
                                                                        repeat_probability = 0.0))
        sticky_seq = generate(sticky, n; rng = StableRNG(322))
        iid_seq = generate(iidlike, n; rng = StableRNG(322))
        sticky_repeats = count(sticky_seq[i] == sticky_seq[i - 1] for i in 2:n) / (n - 1)
        iid_repeats = count(iid_seq[i] == iid_seq[i - 1] for i in 2:n) / (n - 1)
        @test sticky_repeats > iid_repeats + 0.2
    end

    @testset "target_marginal" begin
        g = DyadicLAMP(0.5, [:a, :b], [0.2, 0.8]; d = 10_000)
        @test target_marginal(g) == [0.2, 0.8]
    end
end

@testset "LAMP (MB1)" begin

    @testset "Constructor — valid" begin
        g = LAMP(0.5, [:a, :b, :c])
        @test g.beta == 0.5
        @test g.d    == 1000
        @test g.epsilon == 0.01
        @test isapprox(sum(g.marginal), 1.0)
        @test isapprox(sum(g.weights), 1.0)
        @test g.weights[1] > g.weights[2] > g.weights[end]
        @test g.transition_matrix == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

        g2 = LAMP(0.3, [:x, :y]; d = 200, epsilon = 0.2)
        @test g2.d == 200
        @test g2.epsilon == 0.2

        P = lamp_repeat_transition([0.4, 0.6]; repeat_probability = 0.7)
        g3 = LAMP(0.4, [:x, :y], [0.4, 0.6]; transition_matrix = P)
        @test g3.transition_matrix == P
    end

    @testset "Constructor — argument errors" begin
        @test_throws ArgumentError LAMP(0.0,  [:a, :b])
        @test_throws ArgumentError LAMP(1.0,  [:a, :b])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b], [0.3, 0.5])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b], [-0.1, 1.1])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b], [0.5, Inf])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b, :c], [0.5, 0.5])
        @test_throws ArgumentError LAMP(0.5,  [:a, :b]; d = 0)
        @test_throws ArgumentError LAMP(0.5,  [:a, :b]; epsilon = -0.1)
        @test_throws ArgumentError LAMP(0.5,  [:a, :b]; epsilon = 1.1)
        @test_throws ArgumentError LAMP(0.5,  [:a, :a])
        @test_throws ArgumentError LAMP(0.5, [:a, :b]; transition_matrix = [0.5 0.5])
        @test_throws ArgumentError LAMP(0.5, [:a, :b]; transition_matrix = [0.5 0.6; 0.5 0.4])
        @test_throws ArgumentError lamp_repeat_transition([0.5, 0.5];
                                                          repeat_probability = -0.1)
        @test_throws ArgumentError lamp_repeat_transition([0.5, 0.5];
                                                          repeat_probability = 1.1)
    end

    @testset "lamp_repeat_transition" begin
        P = lamp_repeat_transition([0.25, 0.75]; repeat_probability = 0.8)
        @test P ≈ [0.85 0.15; 0.05 0.95]
        @test all(sum(P; dims = 2) .≈ 1.0)
        @test stationary_distribution(P) ≈ [0.25, 0.75]
    end

    @testset "generate — output type and length" begin
        g   = LAMP(0.5, ['a', 'b', 'c']; d = 100)
        seq = generate(g, 2_000; rng = StableRNG(10))
        @test length(seq) == 2_000
        @test eltype(seq) == Char
        @test all(c ∈ ('a', 'b', 'c') for c in seq)
    end

    @testset "generate — supports d greater than n without random prehistory" begin
        g = LAMP(0.5, [:a, :b], [1.0, 0.0]; d = 10_000, epsilon = 0.0)
        seq = generate(g, 20; rng = StableRNG(123))
        @test seq == fill(:a, 20)
    end

    @testset "generate — transition matrix controls repeated symbols" begin
        n = 2_000
        sticky = LAMP(0.5, [:a, :b], [0.5, 0.5];
                      d = 200,
                      epsilon = 0.0,
                      transition_matrix = lamp_repeat_transition([0.5, 0.5];
                                                                 repeat_probability = 0.95))
        iidlike = LAMP(0.5, [:a, :b], [0.5, 0.5];
                       d = 200,
                       epsilon = 0.0,
                       transition_matrix = lamp_repeat_transition([0.5, 0.5];
                                                                  repeat_probability = 0.0))
        sticky_seq = generate(sticky, n; rng = StableRNG(321))
        iid_seq = generate(iidlike, n; rng = StableRNG(321))
        sticky_repeats = count(sticky_seq[i] == sticky_seq[i - 1] for i in 2:n) / (n - 1)
        iid_repeats = count(iid_seq[i] == iid_seq[i - 1] for i in 2:n) / (n - 1)
        @test sticky_repeats > iid_repeats + 0.2
    end

    @testset "generate — all symbols reachable" begin
        # LAMP with large d exhibits very long mixing times due to LRD —
        # the sample fraction for individual symbols can deviate wildly from the
        # target marginal in finite sequences.  We test only that all symbols are
        # reachable, i.e., at least one occurrence of each symbol appears across
        # multiple independent short runs.
        g = LAMP(0.7, [:a, :b, :c]; d = 50)
        observed = Set{Symbol}()
        for seed in 1:20
            seq = generate(g, 500; rng = StableRNG(seed))
            union!(observed, seq)
        end
        @test :a ∈ observed
        @test :b ∈ observed
        @test :c ∈ observed
    end

    @testset "generate — rejects n < 1" begin
        g = LAMP(0.5, [:a])
        @test_throws ArgumentError generate(g, 0)
    end

    @testset "ACF is positive and decreasing (statistical)" begin
        g   = LAMP(0.3, [:a, :b]; d = 500, epsilon = 0.005)
        seq = generate(g, 20_000; rng = StableRNG(200))
        x   = Float64.(seq .== :a)

        acf1  = _acf(x, 1)
        acf10 = _acf(x, 10)
        acf50 = _acf(x, 50)

        @test acf1  > 0
        @test acf10 > 0
        @test acf1  > acf10 > acf50
    end

    @testset "target_marginal" begin
        g = LAMP(0.5, [:a, :b], [0.2, 0.8]; d = 50)
        @test target_marginal(g) == [0.2, 0.8]
    end

end
