using Test
using Random
using S5

# Shared seeded RNG — each test file that uses `rng` gets a fresh state
# because Julia re-evaluates each include, so we set it at the top level.
const rng = MersenneTwister(42)

@testset "S5.jl" begin
    include("test_utils.jl")
    include("test_pb1.jl")
    include("test_mb1.jl")
    include("test_mb3.jl")
    include("test_io.jl")
end
