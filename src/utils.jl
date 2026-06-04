"""
    quantize_to_symbols(x, alphabet, marginal) -> Vector

Map a real-valued sequence `x` to symbols from `alphabet` by sample-quantile
thresholding, targeting marginal symbol probabilities `marginal`.

Thresholds are placed at the sample quantiles of `x` corresponding to the
cumulative `marginal` probabilities, so the empirical symbol frequencies exactly
match `marginal` for any finite `n`.

# Arguments
- `x::AbstractVector{<:Real}`: real-valued input sequence.
- `alphabet`: ordered collection of `k` symbols.
- `marginal::AbstractVector{<:Real}`: target symbol probabilities (length `k`,
  must sum to 1).

# Examples
```julia
julia> x = randn(1000)
julia> s = quantize_to_symbols(x, [:L, :M, :H], [0.25, 0.5, 0.25])
julia> count(==(:M), s) / 1000 ≈ 0.5
true
```
"""
function quantize_to_symbols(x::AbstractVector{<:Real}, alphabet,
                              marginal::AbstractVector{<:Real})
    k = length(alphabet)
    length(marginal) == k ||
        throw(ArgumentError(
            "marginal length $(length(marginal)) ≠ alphabet length $k"))
    isapprox(sum(marginal), 1.0; atol = 1e-8) ||
        throw(ArgumentError("marginal must sum to 1, got $(sum(marginal))"))
    cum_p      = cumsum(marginal)
    thresholds = quantile(x, @view(cum_p[1:(k - 1)]))
    return [alphabet[searchsortedfirst(thresholds, xi)] for xi in x]
end

"""
    weighted_sample(rng, weights) -> Int

Draw an index from `1:length(weights)` with probability proportional to `weights`.
Uses the sequential inverse-CDF method; O(k) per call.
"""
function weighted_sample(rng::AbstractRNG, weights::AbstractVector{<:Real})
    u    = rand(rng) * sum(weights)
    cumw = 0.0
    for (i, w) in enumerate(weights)
        cumw += w
        u ≤ cumw && return i
    end
    return length(weights)   # numerical safety for floating-point edge cases
end
