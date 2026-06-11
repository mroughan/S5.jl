"""
    WaveletMarkov(H, alphabet, transition_matrices;
                  regime_weights = uniform, cascade_depth = auto)

Property-based LRD symbol-sequence generator (PB3): multiscale cascade driving
a Markov state machine.

`WaveletMarkov` generates a latent multiscale driver using a simple Haar-style
Gaussian cascade. The driver is rank-binned into regimes, and each regime selects
one Markov transition matrix over `alphabet`.

# Arguments
- `H::Real`: Hurst parameter for the latent multiscale driver, `H ∈ (0.5, 1.0)`.
- `alphabet`: ordered collection of unique symbols.
- `transition_matrices`: vector of row-stochastic `k × k` matrices, one per regime.

# Keyword Arguments
- `regime_weights`: target fraction of time spent in each regime. Defaults to
  uniform over regimes.
- `cascade_depth::Int = 0`: number of dyadic cascade levels. `0` means choose
  `floor(log2(n))` at generation time.

# Complexity
O(n log n + n k) time with O(n + R k²) memory.

# Notes
This is a pragmatic PB3 implementation: it uses a Haar-like multiscale driver
rather than a full calibrated wavelet synthesis package. The important interface
property is present: local bigram structure is controlled by explicit Markov
matrices while a multiscale latent process controls regime persistence.

Symbol-level ACF and spectrum diagnostics only see this regime persistence when
the regimes have different observable stationary symbol distributions. If every
regime has the same stationary marginal, the latent multiscale process may be
mostly hidden from one-hot symbol diagnostics.
"""
struct WaveletMarkov{A, W <: AbstractVector{<:Real}} <: LRDGenerator
    H                   :: Float64
    alphabet            :: A
    transition_matrices :: Vector{Matrix{Float64}}
    regime_weights      :: W
    cascade_depth       :: Int

    function WaveletMarkov{A, W}(H::Float64, alphabet::A,
                                 transition_matrices::Vector{Matrix{Float64}},
                                 regime_weights::W,
                                 cascade_depth::Int) where {A, W <: AbstractVector{<:Real}}
        (0.5 < H < 1.0) ||
            throw(ArgumentError("H must be in (0.5, 1.0), got $H"))
        validate_alphabet(alphabet)
        k = length(alphabet)
        isempty(transition_matrices) &&
            throw(ArgumentError("transition_matrices must be non-empty"))
        all(size(P) == (k, k) for P in transition_matrices) ||
            throw(ArgumentError("each transition matrix must have size ($k, $k)"))
        length(regime_weights) == length(transition_matrices) ||
            throw(ArgumentError("regime_weights length must match number of regimes"))
        cascade_depth ≥ 0 ||
            throw(ArgumentError("cascade_depth must be non-negative"))
        new{A, W}(H, alphabet, transition_matrices, regime_weights, cascade_depth)
    end
end

function WaveletMarkov(H::Real, alphabet,
                       transition_matrices::AbstractVector{<:AbstractMatrix{<:Real}};
                       regime_weights::AbstractVector{<:Real} =
                           fill(1.0 / length(transition_matrices),
                                length(transition_matrices)),
                       cascade_depth::Int = 0)
    Ps = [validate_transition_matrix(P, "transition_matrices[$i]")
          for (i, P) in enumerate(transition_matrices)]
    w = validate_probability_vector(regime_weights, "regime_weights")
    WaveletMarkov{typeof(alphabet), typeof(w)}(Float64(H), alphabet, Ps, w,
                                               cascade_depth)
end

"""
    WaveletMarkov(H, specs; regime_weights = uniform, cascade_depth = auto)

Construct a [`WaveletMarkov`](@ref) generator from one [`MarkovSpec`](@ref) per
latent regime. All specifications must use the same ordered alphabet.
"""
function WaveletMarkov(H::Real, specs::AbstractVector{<:MarkovSpec};
                       regime_weights::AbstractVector{<:Real} =
                           fill(1.0 / length(specs), length(specs)),
                       cascade_depth::Int = 0)
    alphabet, Ps = unpack_markov_specs(specs)
    return WaveletMarkov(H, alphabet, Ps; regime_weights, cascade_depth)
end

function Base.show(io::IO, g::WaveletMarkov)
    print(io, "WaveletMarkov{$(typeof(g.alphabet)), $(typeof(g.regime_weights))}",
          "(H=$(g.H), k=$(length(g.alphabet)), R=$(length(g.transition_matrices)))")
end

"""
    generate(g::WaveletMarkov, n; rng) -> Vector

Generate `n` symbols from a [`WaveletMarkov`](@ref) generator.
"""
function generate(g::WaveletMarkov, n::Int; rng::AbstractRNG = Random.default_rng())
    n ≥ 2 || throw(ArgumentError("n must be ≥ 2, got $n"))

    driver = _haar_lrd_driver(n, g.H, g.cascade_depth, rng)
    regimes = quantize_to_symbols(driver, collect(1:length(g.transition_matrices)),
                                  g.regime_weights)
    symbol = weighted_sample(rng, target_marginal(g))
    result = Vector{eltype(g.alphabet)}(undef, n)

    @inbounds for t in 1:n
        P = g.transition_matrices[regimes[t]]
        symbol = weighted_sample(rng, @view P[symbol, :])
        result[t] = g.alphabet[symbol]
    end

    return result
end

function _haar_lrd_driver(n::Int, H::Float64, cascade_depth::Int,
                          rng::AbstractRNG)
    depth = cascade_depth == 0 ? max(1, floor(Int, log2(n))) : cascade_depth
    x = zeros(Float64, n)

    for level in 0:depth
        block = 2^level
        σ = block^(H - 0.5)
        pos = 1
        while pos ≤ n
            val = σ * randn(rng)
            stop = min(n, pos + block - 1)
            @inbounds for i in pos:stop
                x[i] += val
            end
            pos += block
        end
    end

    x .-= mean(x)
    sx = std(x)
    sx > 0 && (x ./= sx)
    return x
end
