@testset "Uniform method factory" begin

    @testset "method metadata" begin
        ids = method_ids()
        @test ids == (:PB1, :PB2, :PB3, :PB4,
                      :MB1a, :MB1b, :MB1c, :MB2, :MB3, :MB4, :MB5)
        @test method_ids(; family = :property_based) == (:PB1, :PB2, :PB3, :PB4)
        @test method_ids(; family = :model_based) == (:MB1a, :MB1b, :MB1c,
                                                       :MB2, :MB3, :MB4, :MB5)
        @test_throws ArgumentError method_ids(; family = :other)

        info = method_info(:MB5)
        @test info.id == :MB5
        @test info.type_name == :DuplicationMutation
        @test info.family == :model_based
        @test :alpha in keys(info.defaults)
        @test method_info("PB1").type_name == :SpectralFGN
        @test method_info(:SpectralFGN).id == :PB1
        @test length(method_info()) == length(ids)
        @test_throws ArgumentError method_info(:NoSuchMethod)
    end

    @testset "construct standard cases" begin
        alphabet = [:a, :b, :c]
        marginal = [0.2, 0.3, 0.5]
        for id in method_ids()
            g = make_generator(id, alphabet; marginal)
            seq = generate(g, 128; rng = StableRNG(100 + findfirst(==(id), method_ids())))
            @test length(seq) == 128
            @test eltype(seq) == Symbol
            @test all(s in alphabet for s in seq)
        end
    end

    @testset "overrides and aliases" begin
        g1 = make_generator(:PB1, [:a, :b]; H = 0.75)
        @test g1 isa SpectralFGN
        @test g1.H == 0.75

        g2 = make_generator("MB1c", [:a, :b]; beta = 0.5, d = 20)
        @test g2 isa CalibratedAdditiveMarkov
        @test g2.beta == 0.5
        @test g2.d == 20

        g3 = make_generator(:LAMP, [:a, :b]; case = :repeat,
                            repeat_probability = 0.8)
        @test g3 isa LAMP
        @test g3.transition_matrix[1, 1] > g3.transition_matrix[1, 2]
    end

    @testset "helpful errors" begin
        @test_throws ArgumentError make_generator(:PB1, [:a, :b]; HH = 0.8)
        @test_throws ArgumentError make_generator(:PB1, [:a, :b]; case = :persistent)
        @test_throws ArgumentError make_generator(:NoSuchMethod, [:a, :b])
    end

end
