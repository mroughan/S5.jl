"""
    LAMP(beta, alphabet [, marginal]; d = 1000, epsilon = 0.01,
         transition_matrix = identity)

Model-based LRD symbol-sequence generator (MB1a): exact finite-history
Linear-Additive Markov Process.

At each step the probability of the next symbol is a convex combination of
transition-matrix rows selected by the most recent `d` history symbols, mixed
with an optional innovation term:

    q(s) = (1 - epsilon) * Σⱼ wⱼ * P[Xₜ₋ⱼ, s] + epsilon * p(s)

with power-law weights `wⱼ ∝ j^{-(1+β)}`, so the autocovariance decays as a
power law with exponent `β` up to the finite history depth, giving nominal
Hurst parameter `H = (2−β)/2`.

# Arguments
- `beta::Real`: ACF decay exponent, `β ∈ (0, 1)`.
- `alphabet`: ordered collection of symbols.
- `marginal::AbstractVector{<:Real}`: stationary marginal (default: uniform).

# Keyword Arguments
- `d::Int = 1000`: history depth. The effective LRD range is bounded by `d`; for
  finite sequences `d` may exceed `n`. Only observed history contributes; missing
  pre-history mass is assigned to the target marginal.
- `epsilon::Real = 0.01`: marginal innovation probability. Larger values improve
  finite-sample marginal control but weaken history dependence.
- `transition_matrix`: row-stochastic transition matrix over `alphabet`. The
  default is identity, so history symbols tend to copy themselves. Use
  [`lamp_repeat_transition`](@ref) for a simple identity/dyad mixture.

# Complexity
O(n·min(d,n)) time, O(d + n) memory.

# References
Kumar, R., Raghu, M., Sarlos, T., & Tomkins, A. (2017). Linear additive Markov
processes. *WWW '17*, 411–419.

Singh, M., Greenberg, C., & Klakow, D. (2016). The custom decay language model
for long range dependencies. *TSD*, 343–351.

# Examples
```julia
julia> g = LAMP(0.5, [:a, :b, :c]; d = 500, epsilon = 0.02)
LAMP{Vector{Symbol}, Vector{Float64}}(β=0.5, k=3, d=500, ε=0.02)

julia> seq = generate(g, 5000; rng = MersenneTwister(42))
julia> length(seq) == 5000 && eltype(seq) == Symbol
true
```
"""
struct LAMP{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    beta              :: Float64
    alphabet          :: A
    marginal          :: M
    d                 :: Int
    weights           :: Vector{Float64}
    epsilon           :: Float64
    transition_matrix :: Matrix{Float64}

    function LAMP{A, M}(beta::Float64, alphabet::A, marginal::M,
                         d::Int, weights::Vector{Float64},
                         epsilon::Float64,
                         transition_matrix::Matrix{Float64}) where {A, M <: AbstractVector{<:Real}}
        (0.0 < beta < 1.0) ||
            throw(ArgumentError("beta must be in (0, 1), got $beta"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        length(marginal) == k ||
            throw(ArgumentError(
                "marginal length $(length(marginal)) ≠ alphabet length $k"))
        d ≥ 1 || throw(ArgumentError("d must be ≥ 1, got $d"))
        0.0 ≤ epsilon ≤ 1.0 ||
            throw(ArgumentError("epsilon must be in [0, 1], got $epsilon"))
        size(transition_matrix) == (k, k) ||
            throw(ArgumentError(
                "transition_matrix must have size ($k, $k), got $(size(transition_matrix))"))
        new{A, M}(beta, alphabet, marginal, d, weights, epsilon, transition_matrix)
    end
end

function LAMP(beta::Real, alphabet,
              marginal::AbstractVector{<:Real} =
                  fill(1.0 / length(alphabet), length(alphabet));
              d::Int = 1000,
              epsilon::Real = 0.01,
              transition_matrix::Union{Nothing, AbstractMatrix{<:Real}} = nothing)
    m  = validate_probability_vector(marginal, "marginal")
    w  = [j^(-(1.0 + Float64(beta))) for j in 1:d]
    w ./= sum(w)
    P = transition_matrix === nothing ?
        Matrix{Float64}(I, length(alphabet), length(alphabet)) :
        validate_transition_matrix(transition_matrix, "transition_matrix")
    LAMP{typeof(alphabet), typeof(m)}(Float64(beta), alphabet, m, d, w,
                                      Float64(epsilon), P)
end

function Base.show(io::IO, g::LAMP)
    print(io, "LAMP{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(β=$(g.beta), k=$(length(g.alphabet)), d=$(g.d), ε=$(g.epsilon))")
end

"""
    lamp_repeat_transition(marginal; repeat_probability = 0.8) -> Matrix{Float64}

Construct a row-stochastic identity/dyad transition matrix for [`LAMP`](@ref).

The matrix is

    P[i, j] = repeat_probability * 𝟏[i = j] +
              (1 - repeat_probability) * marginal[j]

so larger `repeat_probability` makes the process more likely to repeat the
state selected from history, while the dyad term pulls rows back toward the
target marginal.

# Examples
```julia
julia> P = lamp_repeat_transition([0.25, 0.75]; repeat_probability = 0.8)
2×2 Matrix{Float64}:
 0.85  0.15
 0.05  0.95
```
"""
function lamp_repeat_transition(marginal::AbstractVector{<:Real};
                                repeat_probability::Real = 0.8)
    p = validate_probability_vector(marginal, "marginal")
    ρ = Float64(repeat_probability)
    0.0 ≤ ρ ≤ 1.0 ||
        throw(ArgumentError("repeat_probability must be in [0, 1], got $repeat_probability"))
    k = length(p)
    P = repeat(reshape((1 - ρ) .* p, 1, k), k, 1)
    @inbounds for i in 1:k
        P[i, i] += ρ
    end
    return P
end

"""
    generate(g::LAMP, n; rng) -> Vector

Generate `n` symbols from a [`LAMP`](@ref) generator using observed history.

The next symbol is drawn from a probability vector formed by weighting transition
rows selected by previous symbols. If `d` exceeds the available observed history,
the missing pre-history weight is assigned to `g.marginal`. This makes `d > n`
well-defined without inventing an unobserved random history.
"""
function generate(g::LAMP, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    k = length(g.alphabet)
    result_indices = Vector{Int}(undef, n)
    result = Vector{eltype(g.alphabet)}(undef, n)
    q = Vector{Float64}(undef, k)

    @inbounds for t in 1:n
        fill!(q, 0.0)
        maxhist = min(g.d, t - 1)
        observed_weight = 0.0
        for j in 1:maxhist
            wj = g.weights[j]
            observed_weight += wj
            row = result_indices[t - j]
            for s in 1:k
                q[s] += wj * g.transition_matrix[row, s]
            end
        end
        missing_weight = 1 - observed_weight
        if missing_weight > 0
            for s in 1:k
                q[s] += missing_weight * g.marginal[s]
            end
        end
        if g.epsilon > 0
            for s in 1:k
                q[s] = (1 - g.epsilon) * q[s] + g.epsilon * g.marginal[s]
            end
        end

        idx = weighted_sample(rng, q)
        result_indices[t] = idx
        result[t] = g.alphabet[idx]
    end

    return result
end

"""
    DyadicLAMP(beta, alphabet [, marginal]; d = 1_000_000, epsilon = 0.01,
               transition_matrix = identity)

Scalable dyadic-bucket approximation to [`LAMP`](@ref) (MB1b).

`DyadicLAMP` uses the same power-law lag weights and transition-matrix control
as exact LAMP, but compresses observed history into age buckets
`1`, `2:3`, `4:7`, and so on. Each bucket contributes its total power-law weight
times the empirical symbol mix in that bucket. Missing pre-history mass is
assigned to the target marginal.

This is a finite-sequence approximation for large `d` and long sequences; it is
not an exact replacement for [`LAMP`](@ref).

# Complexity
O(n · k · log(n) · log(min(d, n))) time and O(n · k) memory, where `k` is the
alphabet size.

# Examples
```julia
julia> p = [0.4, 0.6]
julia> P = lamp_repeat_transition(p; repeat_probability = 0.8)
julia> g = DyadicLAMP(0.5, [:a, :b], p; d = 10_000, transition_matrix = P)
DyadicLAMP{Vector{Symbol}, Vector{Float64}}(β=0.5, k=2, d=10000, ε=0.01)

julia> length(generate(g, 1000; rng = MersenneTwister(1))) == 1000
true
```
"""
struct DyadicLAMP{A, M <: AbstractVector{<:Real}} <: LRDGenerator
    beta              :: Float64
    alphabet          :: A
    marginal          :: M
    d                 :: Int
    weights           :: Vector{Float64}
    prefix_weights    :: Vector{Float64}
    epsilon           :: Float64
    transition_matrix :: Matrix{Float64}

    function DyadicLAMP{A, M}(beta::Float64, alphabet::A, marginal::M,
                              d::Int, weights::Vector{Float64},
                              prefix_weights::Vector{Float64},
                              epsilon::Float64,
                              transition_matrix::Matrix{Float64}) where {A, M <: AbstractVector{<:Real}}
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
        length(prefix_weights) == d + 1 ||
            throw(ArgumentError("prefix_weights length must equal d + 1"))
        0.0 ≤ epsilon ≤ 1.0 ||
            throw(ArgumentError("epsilon must be in [0, 1], got $epsilon"))
        size(transition_matrix) == (k, k) ||
            throw(ArgumentError(
                "transition_matrix must have size ($k, $k), got $(size(transition_matrix))"))
        new{A, M}(beta, alphabet, marginal, d, weights, prefix_weights,
                  epsilon, transition_matrix)
    end
end

function _lamp_weights(beta::Real, d::Int)
    d ≥ 1 || throw(ArgumentError("d must be ≥ 1, got $d"))
    w = [j^(-(1.0 + Float64(beta))) for j in 1:d]
    w ./= sum(w)
    return w
end

function _lamp_prefix_weights(weights::AbstractVector{<:Real})
    prefix = zeros(Float64, length(weights) + 1)
    @inbounds for i in eachindex(weights)
        prefix[i + 1] = prefix[i] + weights[i]
    end
    prefix[end] = 1.0
    return prefix
end

function _lamp_transition_matrix(alphabet, transition_matrix)
    return transition_matrix === nothing ?
        Matrix{Float64}(I, length(alphabet), length(alphabet)) :
        validate_transition_matrix(transition_matrix, "transition_matrix")
end

function DyadicLAMP(beta::Real, alphabet,
                    marginal::AbstractVector{<:Real} =
                        fill(1.0 / length(alphabet), length(alphabet));
                    d::Int = 1_000_000,
                    epsilon::Real = 0.01,
                    transition_matrix::Union{Nothing, AbstractMatrix{<:Real}} = nothing)
    m = validate_probability_vector(marginal, "marginal")
    w = _lamp_weights(beta, d)
    prefix = _lamp_prefix_weights(w)
    P = _lamp_transition_matrix(alphabet, transition_matrix)
    return DyadicLAMP{typeof(alphabet), typeof(m)}(
        Float64(beta), alphabet, m, d, w, prefix, Float64(epsilon), P)
end

function Base.show(io::IO, g::DyadicLAMP)
    print(io, "DyadicLAMP{$(typeof(g.alphabet)), $(typeof(g.marginal))}",
          "(β=$(g.beta), k=$(length(g.alphabet)), d=$(g.d), ε=$(g.epsilon))")
end

function _fenwick_add!(tree::Matrix{Int}, symbol::Int, index::Int)
    n = size(tree, 2)
    while index ≤ n
        tree[symbol, index] += 1
        index += index & -index
    end
    return tree
end

function _fenwick_prefix(tree::Matrix{Int}, symbol::Int, index::Int)
    total = 0
    while index > 0
        total += tree[symbol, index]
        index -= index & -index
    end
    return total
end

function _fenwick_range_counts!(counts::Vector{Int}, tree::Matrix{Int},
                                lo::Int, hi::Int)
    fill!(counts, 0)
    hi < lo && return counts
    @inbounds for symbol in eachindex(counts)
        counts[symbol] = _fenwick_prefix(tree, symbol, hi) -
                         _fenwick_prefix(tree, symbol, lo - 1)
    end
    return counts
end

"""
    generate(g::DyadicLAMP, n; rng) -> Vector

Generate `n` symbols from the dyadic-bucket LAMP approximation.

The previous observations are summarized in dyadic age buckets. Bucket
`start:stop` contributes the total power-law mass for those lags multiplied by
the empirical symbol distribution in that bucket.
"""
function generate(g::DyadicLAMP, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 1 || throw(ArgumentError("n must be ≥ 1, got $n"))
    k = length(g.alphabet)
    tree = zeros(Int, k, n)
    counts = zeros(Int, k)
    q = Vector{Float64}(undef, k)
    result_indices = Vector{Int}(undef, n)
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        fill!(q, 0.0)
        maxhist = min(g.d, t - 1)
        observed_weight = 0.0
        start = 1
        while start ≤ maxhist
            stop = min(2 * start - 1, maxhist)
            bucket_weight = g.prefix_weights[stop + 1] - g.prefix_weights[start]
            observed_weight += bucket_weight
            lo = t - stop
            hi = t - start
            bucket_len = hi - lo + 1
            _fenwick_range_counts!(counts, tree, lo, hi)
            for row in 1:k
                counts[row] == 0 && continue
                coeff = bucket_weight * counts[row] / bucket_len
                for s in 1:k
                    q[s] += coeff * g.transition_matrix[row, s]
                end
            end
            start *= 2
        end

        missing_weight = 1 - observed_weight
        if missing_weight > 0
            for s in 1:k
                q[s] += missing_weight * g.marginal[s]
            end
        end
        if g.epsilon > 0
            for s in 1:k
                q[s] = (1 - g.epsilon) * q[s] + g.epsilon * g.marginal[s]
            end
        end

        idx = weighted_sample(rng, q)
        result_indices[t] = idx
        result[t] = g.alphabet[idx]
        _fenwick_add!(tree, idx, t)
    end

    return result
end
