"""
    CalibratedAdditiveMarkov(beta, alphabet [, marginal]; d = 1000,
                             strength = 0.8)

Model-based symbolic generator (MB1c) using a centered additive Markov memory
function.

At each step, the next-symbol probabilities are

    q(s) = p(s) + strength * Σⱼ wⱼ * (1[X[t-j] = s] - p(s)),

where `p` is the target marginal and `wⱼ ∝ j^(-beta)` over `1:d`. The centered
terms sum to zero across symbols, so `q` remains a probability vector for
`strength ∈ [0, 1]`. Larger `strength` gives more weight to observed history;
`strength = 0` gives iid draws from `marginal`.

This generator is related to additive Markov-chain memory-function models. It is
finite-history and does not claim exact asymptotic LRD beyond the configured
memory depth `d`.

# Arguments
- `beta::Real`: nominal memory-function decay exponent, `beta ∈ (0, 1)`.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: target marginal probabilities (default:
  uniform).

# Keyword Arguments
- `d::Int = 1000`: history depth and finite memory cutoff.
- `strength::Real = 0.8`: history coupling strength in `[0, 1]`.

# Complexity
O(n·min(d,n)) time, O(n) memory.

# References
Melnyk, S. S., Usatenko, O. V., & Yampol'skii, V. A. (2006). Memory functions
of the additive Markov chains: applications to complex dynamic systems.
*Physica A* 361, 405-415.

Mayzelis, Z. A., Apostolov, S. S., Melnyk, S. S., Usatenko, O. V., &
Yampol'skii, V. A. (2006). Additive N-step Markov chains as prototype model of
symbolic stochastic dynamical systems with long-range correlations.

# Examples
```julia
julia> g = CalibratedAdditiveMarkov(0.4, [:x, :y], [0.3, 0.7]; d = 200)
CalibratedAdditiveMarkov{Vector{Symbol}, Vector{Float64}}(β=0.4, k=2, d=200, strength=0.8)

julia> length(generate(g, 1000; rng = MersenneTwister(1))) == 1000
true
```
"""
struct CalibratedAdditiveMarkov{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    beta     :: Float64
    alphabet :: A
    marginal :: M
    d        :: Int
    weights  :: Vector{Float64}
    strength :: Float64

    function CalibratedAdditiveMarkov{A, M}(beta::Float64, alphabet::A,
                                            marginal::M, d::Int,
                                            weights::Vector{Float64},
                                            strength::Float64) where {A, M <: AbstractVector{<:Real}}
        (0.0 < beta < 1.0) ||
            throw(ArgumentError("beta must be in (0, 1), got $beta"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        d ≥ 1 || throw(ArgumentError("d must be ≥ 1, got $d"))
        length(weights) == d ||
            throw(ArgumentError("weights length must equal d"))
        0.0 ≤ strength ≤ 1.0 ||
            throw(ArgumentError("strength must be in [0, 1], got $strength"))
        new{A, M}(beta, alphabet, marginal, d, weights, strength)
    end
end

function CalibratedAdditiveMarkov(beta::Real, alphabet,
                                  marginal::AbstractVector{<:Real} =
                                      fill(1.0 / length(alphabet), length(alphabet));
                                  d::Int = 1000,
                                  strength::Real = 0.8)
    m = validate_probability_vector(marginal, "marginal")
    d ≥ 1 || throw(ArgumentError("d must be ≥ 1, got $d"))
    w = [j^(-Float64(beta)) for j in 1:d]
    w ./= sum(w)
    CalibratedAdditiveMarkov{typeof(alphabet), typeof(m)}(
        Float64(beta), alphabet, m, d, w, Float64(strength))
end

function Base.show(io::IO, g::CalibratedAdditiveMarkov)
    print(io, "CalibratedAdditiveMarkov{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(β=$(g.beta), k=$(length(g.alphabet)), d=$(g.d), strength=$(g.strength))")
end

"""
    generate(g::CalibratedAdditiveMarkov, n; rng) -> Vector

Generate `n` symbols from a [`CalibratedAdditiveMarkov`](@ref) generator.
"""
function generate(g::CalibratedAdditiveMarkov, n::Int;
                  rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    k = length(g.alphabet)
    result_indices = Vector{Int}(undef, n)
    result = Vector{eltype(g.alphabet)}(undef, n)
    q = Vector{Float64}(undef, k)

    @inbounds for t in 1:n
        q .= (1 - g.strength) .* g.marginal
        maxhist = min(g.d, t - 1)
        observed_weight = 0.0
        for j in 1:maxhist
            wj = g.weights[j]
            observed_weight += wj
            idx = result_indices[t - j]
            for s in 1:k
                q[s] += g.strength * wj * (s == idx ? 1.0 : 0.0)
            end
        end
        missing_weight = 1 - observed_weight
        if missing_weight > 0
            for s in 1:k
                q[s] += g.strength * missing_weight * g.marginal[s]
            end
        end

        idx = weighted_sample(rng, q)
        result_indices[t] = idx
        result[t] = g.alphabet[idx]
    end

    return result
end

"""
    DuplicationMutation(alpha, alphabet [, marginal]; mutation_probability = 0.01,
                        seed_length = 64, max_block_length = 4096)

Model-based symbolic growth generator (MB5) using copy-and-mutate block
duplication.

Generation starts with `seed_length` iid symbols from `marginal`. The sequence
then grows one symbol at a time by choosing a power-law copy distance, copying
the symbol at that lag, and mutating the copied symbol independently with
probability `mutation_probability`. Copy distances are drawn from a truncated
power law `P(D = ell) ∝ ell^(-alpha)` over the available history, capped by
`max_block_length`.

This is a finite-sequence symbolic analogue of expansion-modification and
duplication-mutation ideas. It is naturally DNA-like, but it does not provide
direct bigram control or an exact Hurst-parameter guarantee. The power-law copy
distance is the part that gives a pathway to broad lag dependence. Earlier
block-copy variants with uniformly chosen source blocks mostly created local
duplication patches rather than a power-law autocorrelation curve.

# Arguments
- `alpha::Real`: copy-distance exponent, `alpha > 1`; use values near `1`
  for slower empirical decay and broader dependence.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: mutation replacement and seed marginal
  probabilities (default: uniform).

# Keyword Arguments
- `mutation_probability::Real = 0.01`: per-symbol mutation probability.
- `seed_length::Int = 64`: iid prefix length before copy-mutate growth.
- `max_block_length::Int = 4096`: legacy keyword naming the maximum copy
  distance/backward memory window.

# Complexity
O(n log d + d) time and O(n + d) memory, where `d = max_block_length`.

# References
Li, W. (1991). Expansion-modification systems: a model for spatial 1/f spectra.
*Physical Review A* 43, 5240-5260.

Li, W., Marr, T. G., & Kaneko, K. (1994). Understanding long-range correlations
in DNA sequences. *Physica D* 75, 392-416.

# Examples
```julia
julia> g = DuplicationMutation(1.5, ['A', 'C', 'G', 'T']; mutation_probability = 0.02)
DuplicationMutation{Vector{Char}, Vector{Float64}}(α=1.5, k=4, μ=0.02, seed=64, max_block=4096)

julia> length(generate(g, 500; rng = MersenneTwister(1))) == 500
true
```
"""
struct DuplicationMutation{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    alpha                :: Float64
    alphabet             :: A
    marginal             :: M
    mutation_probability :: Float64
    seed_length          :: Int
    max_block_length     :: Int

    function DuplicationMutation{A, M}(alpha::Float64, alphabet::A,
                                       marginal::M,
                                       mutation_probability::Float64,
                                       seed_length::Int,
                                       max_block_length::Int) where {A, M <: AbstractVector{<:Real}}
        alpha > 1.0 || throw(ArgumentError("alpha must be > 1, got $alpha"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        0.0 ≤ mutation_probability ≤ 1.0 ||
            throw(ArgumentError(
                "mutation_probability must be in [0, 1], got $mutation_probability"))
        seed_length ≥ 1 ||
            throw(ArgumentError("seed_length must be ≥ 1, got $seed_length"))
        max_block_length ≥ 1 ||
            throw(ArgumentError(
                "max_block_length must be ≥ 1, got $max_block_length"))
        new{A, M}(alpha, alphabet, marginal, mutation_probability, seed_length,
                  max_block_length)
    end
end

function DuplicationMutation(alpha::Real, alphabet,
                             marginal::AbstractVector{<:Real} =
                                 fill(1.0 / length(alphabet), length(alphabet));
                             mutation_probability::Real = 0.01,
                             seed_length::Int = 64,
                             max_block_length::Int = 4096)
    m = validate_probability_vector(marginal, "marginal")
    DuplicationMutation{typeof(alphabet), typeof(m)}(
        Float64(alpha), alphabet, m, Float64(mutation_probability),
        seed_length, max_block_length)
end

function Base.show(io::IO, g::DuplicationMutation)
    print(io, "DuplicationMutation{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(α=$(g.alpha), k=$(length(g.alphabet)), μ=$(g.mutation_probability), ",
          "seed=$(g.seed_length), max_block=$(g.max_block_length))")
end

function _power_law_cdf(exponent::Float64, maxlen::Int)
    cdf = Vector{Float64}(undef, maxlen)
    total = 0.0
    @inbounds for j in 1:maxlen
        total += j^(-exponent)
        cdf[j] = total
    end
    return cdf
end

function _sample_power_law_distance(rng::AbstractRNG, cdf::Vector{Float64},
                                    maxdistance::Int)
    u = rand(rng) * cdf[maxdistance]
    return searchsortedfirst(cdf, u, 1, maxdistance, Base.Order.Forward)
end

"""
    generate(g::DuplicationMutation, n; rng) -> Vector

Generate `n` symbols from a [`DuplicationMutation`](@ref) generator.
"""
function generate(g::DuplicationMutation, n::Int;
                  rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    result_indices = Vector{Int}(undef, n)
    result = Vector{eltype(g.alphabet)}(undef, n)
    copy_cdf = _power_law_cdf(g.alpha, min(g.max_block_length, max(1, n - 1)))

    seed_n = min(n, g.seed_length)
    @inbounds for t in 1:seed_n
        idx = weighted_sample(rng, g.marginal)
        result_indices[t] = idx
        result[t] = g.alphabet[idx]
    end

    pos = seed_n + 1
    @inbounds while pos ≤ n
        maxdistance = min(length(copy_cdf), pos - 1)
        distance = _sample_power_law_distance(rng, copy_cdf, maxdistance)
        idx = result_indices[pos - distance]
        if rand(rng) < g.mutation_probability
            idx = weighted_sample(rng, g.marginal)
        end
        result_indices[pos] = idx
        result[pos] = g.alphabet[idx]
        pos += 1
    end

    return result
end
