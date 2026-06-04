# S5.jl — Implementation TODO

This file tracks the planned implementation of all six LRD symbol-sequence synthesis
methods. Methods are grouped into property-based and model-based categories as described
in the ARC Discovery Grant proposal (Roughan & Willinger 2023).

For each method the following must be completed:
- [ ] Core generator function(s) with docstrings and examples
- [ ] Unit tests validating basic output (type, length, alphabet membership)
- [ ] Statistical validation test: estimate $H$ from generated sequence and compare to target
- [ ] Benchmark: measure wall time and memory for $n \in \{10^4, 10^5, 10^6\}$
- [ ] Update CHANGELOG.md on completion

---

## Priority 1 (implement first)

### PB1 — Spectral fGn + Quantization

**Status:** Not started

**Description.**
Generate fractional Gaussian noise (fGn) with Hurst parameter $H \in (1/2, 1)$ using
the fast spectral (FFT) method, then map the real-valued output to symbols by amplitude
quantization.

**Algorithm.**
1. Construct the fGn power spectrum $S(f) \propto |f|^{-(2H-1)}$ on the DFT grid.
2. Fill with complex Gaussian noise scaled by $\sqrt{S(f)}$, respecting conjugate symmetry.
3. Take the real part of the IFFT to obtain fGn of length $n$.
4. Map to symbols: sort thresholds to match a target marginal distribution $\mathbf{p}$ over
   the alphabet; assign each sample to the bin it falls in.

**Key parameters.**
- `H::Float64` — Hurst parameter, $H \in (1/2, 1)$
- `n::Int` — sequence length
- `alphabet` — ordered collection of symbols (e.g., `['a','b','c']`)
- `marginal::Vector{Float64}` — target symbol probabilities (defaults to uniform)

**References.** Paxson (1997) CCR 27; Dieker (2004) PhD thesis U. Twente.

**Validation.**
- Estimate $H$ from the generated sequence using a spectral log-log regression (the
  periodogram slope at low frequencies) and check it matches the target within tolerance.
- Check marginal symbol frequencies converge to `marginal` as $n \to \infty$.

**Known limitations.**
- Short-range structure (bigrams, etc.) is determined entirely by the quantization grid
  and cannot be prescribed independently.
- The Paxson method is approximate at high frequencies; the circulant embedding
  alternative is exact but requires $n$ to be a power of 2 or uses zero-padding.

---

### MB1 — Linear-Additive Markov Process (LAMP)

**Status:** Not started

**Description.**
Generate a finite-state symbol sequence whose transition probabilities are a weighted
linear combination of the recent history, with power-law weights that enforce a
prescribed ACF decay exponent $\beta$.

**Algorithm.**
1. Choose history depth $d$ (truncation of the infinite sum).
2. Construct weights $w_k = k^{-(1+\beta)} / \sum_{j=1}^{d} j^{-(1+\beta)}$ for $k = 1, \ldots, d$.
3. Initialise with $d$ symbols drawn from the marginal distribution $\mathbf{p}$.
4. At each step $t$, compute the (unnormalised) probability vector
   $\mathbf{q} = \sum_{k=1}^{d} w_k \, \mathbf{e}_{X_{t-k}}$,
   normalise to obtain $\mathbf{q} / \|\mathbf{q}\|_1$, and draw the next symbol.

**Key parameters.**
- `beta::Float64` — ACF decay exponent, $\beta \in (0, 1)$; maps to $H = (2-\beta)/2$
- `n::Int` — sequence length
- `alphabet` — symbol set
- `d::Int` — history depth (default: auto-selected as $\lceil n^{1/(1+\beta)} \rceil$)
- `marginal::Vector{Float64}` — stationary marginal (defaults to uniform)

**References.** Kumar, Raghu, Sarlos & Tomkins (WWW 2017); Singh, Greenberg & Klakow
(TSD 2016) for the CDLM variant.

**Validation.**
- Estimate $\beta$ from the empirical ACF on a log-log plot and compare to target.
- Verify that marginal symbol frequencies match `marginal`.
- Check that the "missing scales" issue does not arise: confirm LRD is visible across
  at least three decades of lag.

**Known limitations.**
- For large alphabets the $O(|\Sigma| \cdot d)$ cost per step can be expensive; a
  sparse or low-rank representation of the weight tensor may be needed.
- The effective $\beta$ depends on the truncation depth $d$; sequences shorter than
  $d^{1+\beta}$ may not show asymptotic LRD.

---

### MB3 — Fractal Symbol Sequence (FSS) via FRP/FSNP

**Status:** Not started

**Description.**
Each symbol in the alphabet is assigned an independent fractal point process governing
the times at which that symbol is emitted. The sequence is the merge of all symbol
streams, with the symbol whose next event time is smallest being the output at each step.

**Algorithm.**
1. For each symbol $s_i$, draw inter-arrival times $\tau \sim \text{Pareto}(\alpha, x_\text{min})$
   with $\alpha \in (1, 2)$ (giving $H = (3-\alpha)/2$). Use a Fractal Shot Noise
   Process (FSNP) rather than a naive FRP to avoid the "missing scales" problem [47].
2. Maintain a priority queue of (next\_event\_time, symbol) pairs.
3. Pop the minimum, emit that symbol, draw the next inter-arrival for that symbol, and
   push back.
4. Repeat for $n$ emissions.

**Key parameters.**
- `alpha::Float64` — Pareto tail index, $\alpha \in (1, 2)$; gives $H = (3-\alpha)/2$
- `n::Int` — number of symbols to generate
- `alphabet` — symbol set
- `rates::Vector{Float64}` — base arrival rates per symbol (controls marginal frequencies)
- `x_min::Float64` — Pareto scale parameter

**References.** Lowen & Teich (Fractals 1995); Ryu & Lowen (Stochastic Models 1998);
Roughan, Yates & Veitch (1999) on the missing-scales pitfall.

**Validation.**
- Confirm that the count-process variance grows faster than linearly (Definition 3 from
  the proposal), and fit the growth exponent.
- Verify the LRD scale range covers at least three decades.
- Check marginal symbol frequencies match `rates / sum(rates)`.

**Known limitations.**
- Independent streams for each symbol make it hard to prescribe joint symbol statistics
  (bigrams); this is a fundamental limitation of the FSS approach.
- Care is required to ensure the FSNP scale range is not smaller than the sequence
  length (the Roughan–Yates–Veitch pitfall).

---

## Priority 2 (implement after Priority 1 and test infrastructure is in place)

### PB2 — Latent Gaussian Categorical Model (LGCM)

**Status:** Not started

**Description.**
A vector of $k$ correlated Gaussian processes, one per symbol, shares an fGn covariance
structure. At each time step the symbol is the argmax of the latent vector. Per-symbol
means shift marginal probabilities.

**Algorithm.**
1. Construct the fGn covariance matrix $\Sigma_{ij} = \tfrac{1}{2}(|i|^{2H} + |j|^{2H} - |i-j|^{2H})$
   for indices $i, j = 1, \ldots, n$ (or use a Toeplitz approximation for large $n$).
2. For each symbol $s_m$, draw a length-$n$ Gaussian vector $\mathbf{z}^{(m)} \sim \mathcal{N}(\mu_m \mathbf{1}, \Sigma)$.
3. At each step $t$, emit symbol $\arg\max_m z^{(m)}_t$.
4. Set means $\mu_m$ to achieve a target marginal via a root-finding step (bisection on
   the resulting multinomial probabilities, using Monte Carlo estimates if needed).

**Key parameters.**
- `H::Float64` — Hurst parameter
- `n::Int` — sequence length
- `alphabet` — symbol set
- `marginal::Vector{Float64}` — target symbol probabilities
- `approx::Symbol` — `:exact` (Cholesky, $O(n^2)$ memory) or `:circulant` (FFT, approximate)

**References.** Gal, Chen & Ghahramani (ICML 2015).

**Validation.**
- Estimate $H$ from the indicator series $\mathbf{1}[X_t = s]$ for each symbol and
  check consistency with input $H$.
- Verify marginal symbol frequencies match `marginal`.

**Known limitations.**
- Exact Cholesky factorisation is $O(n^3)$ time and $O(n^2)$ memory; only feasible for
  short sequences. Circulant/FFT approximation is needed for $n > 10^4$.

---

### PB3 — Wavelet-Cascade Driving a Markov State Machine

**Status:** Not started

**Description.**
A latent LRD intensity signal is synthesised via a wavelet cascade. This signal
continuously selects among a set of Markov transition matrices, so local (bigram)
structure is prescribed by the matrices while LRD is injected at all scales by the
wavelet layer.

**Algorithm.**
1. Generate a wavelet-synthesised LRD signal $\lambda_t$ (using a log-normal or
   multiplicative cascade on wavelet coefficients, as in Roughan, Veitch & Abry 2001).
2. Define $R$ "regimes" with corresponding Markov transition matrices
   $\{P^{(1)}, \ldots, P^{(R)}\}$. Partition the range of $\lambda_t$ into $R$ bins.
3. At each step $t$, look up the regime $r_t = \text{bin}(\lambda_t)$ and draw the next
   symbol from the row $P^{(r_t)}_{X_{t-1}, \cdot}$.

**Key parameters.**
- `H::Float64` — Hurst parameter for the wavelet layer
- `n::Int` — sequence length
- `transition_matrices::Vector{Matrix{Float64}}` — one $|\Sigma| \times |\Sigma|$ matrix per regime
- `n_regimes::Int` — number of regimes $R$ (default: 2)

**References.** Roughan, Veitch & Abry (IEEE/ACM ToN 2001).

**Validation.**
- Estimate $H$ from the generated sequence and compare to the wavelet-layer target.
- Verify that empirical bigram frequencies match those implied by the transition matrices
  and the stationary distribution over regimes.

**Known limitations.**
- The coupling between the wavelet layer and the Markov layer adds parameters that must
  be calibrated; the effective $H$ of the symbol sequence is not identical to the $H$
  of the wavelet signal and must be verified empirically.

---

### MB2 — Heavy-Tailed On/Off Doubly-Stochastic Markov Chain

**Status:** Not started

**Description.**
A Markov chain alternates between two (or more) regimes. Sojourn times follow a Pareto
distribution with tail index $\alpha \in (1, 2)$, giving $H = (3-\alpha)/2$. Within
each regime a distinct (SRD) Markov chain governs symbol emissions.

**Algorithm.**
1. Initialise in regime $r \in \{1, \ldots, R\}$ with a sojourn length
   $L \sim \text{Pareto}(\alpha, L_\text{min})$ rounded to an integer.
2. For $L$ steps, draw symbols from the transition matrix $P^{(r)}$.
3. Draw the next regime $r'$ from a regime-switching matrix $Q$ and a new sojourn $L'$.
4. Repeat until $n$ symbols are generated.

**Key parameters.**
- `alpha::Float64` — Pareto tail index, $\alpha \in (1, 2)$
- `n::Int` — sequence length
- `transition_matrices::Vector{Matrix{Float64}}` — one per regime
- `switching_matrix::Matrix{Float64}` — regime-to-regime transition probabilities $Q$
- `L_min::Float64` — minimum sojourn length

**References.** Garrett & Willinger (ACM Sigcomm 1994); Ryu & Lowen (Stochastic Models 1998).

**Validation.**
- Estimate $H$ from the count-process variance growth and compare to
  $(3-\alpha)/2$.
- Verify that empirical bigram frequencies in each regime match the corresponding
  transition matrix.
- Check for the missing-scales problem: confirm LRD is visible across at least three
  decades of lag.

**Known limitations.**
- The LRD scale range is bounded by the longest sojourn times generated; for short
  sequences this can limit the apparent LRD to fewer decades than desired.

---

## Cross-Cutting Tasks

- [ ] Define a common Julia interface (`generate(method, n, alphabet; kwargs...)`) that
      all six methods implement.
- [ ] Write a shared statistical validation suite (`test/validate_lrd.jl`) that applies
      a spectral $H$ estimator and a count-variance test to any generated sequence.
- [ ] Write benchmarks (`benchmark/benchmarks.jl`) using BenchmarkTools.jl for all methods.
- [ ] Set up Project.toml and Manifest.toml with required dependencies.
- [ ] Create CHANGELOG.md and update it as each method is completed.
- [ ] Add AI disclosure to README.md (done).
