"""
    LRDGenerator

Abstract supertype for all LRD symbol-sequence generators in S5.jl.

Concrete subtypes must implement:

    generate(g::MyGenerator, n::Int; rng::AbstractRNG = Random.default_rng()) -> Vector
"""
abstract type LRDGenerator end

"""
    generate(g, n; rng = Random.default_rng()) -> Vector

Generate a sequence of `n` symbols using LRD generator `g`.

Returns a `Vector` whose element type matches the alphabet of `g`.

# Arguments
- `g::LRDGenerator`: configured generator instance.
- `n::Int`: number of symbols to emit.

# Keyword Arguments
- `rng::AbstractRNG`: random number generator (default: `Random.default_rng()`).

# Examples
```julia
julia> g = SpectralFGN(0.8, ['a', 'b', 'c'])
julia> seq = generate(g, 1000)
julia> length(seq) == 1000
true
```
"""
function generate end
