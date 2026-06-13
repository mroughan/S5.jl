# S5.jl

**Self-Similar Symbols Sequence Synthesis**

S5.jl generates Long-Range Dependent (LRD) sequences of categorical
(non-numerical) symbols. It produces controllable synthetic data for LRD
estimator tests, information-theoretic experiments, and LLM-style neural
sequence models trained on challenging non-language symbolic streams.

## Quick start

```julia
using S5, StableRNGs

rng = StableRNG(42)

# Uniform standard-case factory
g0 = make_generator(:PB1, [:a, :b, :c]; H = 0.8, marginal = [0.2, 0.3, 0.5])
seq0 = generate(g0, 10_000; rng)

method_ids()
method_info(:PB1).defaults

# Property-based: spectral fGn + quantization
g1 = SpectralFGN(0.8, [:a, :b, :c], [0.2, 0.3, 0.5])
seq1 = generate(g1, 10_000; rng)

# Property-based: latent Gaussian categorical model
g2 = LGCM(0.8, [:a, :b, :c], [0.2, 0.3, 0.5])
seq2 = generate(g2, 10_000; rng)

# Property-based: multiscale driver + Markov regimes
P1 = [0.9 0.1; 0.2 0.8]
P2 = [0.3 0.7; 0.6 0.4]
g3 = WaveletMarkov(0.8, [:a, :b], [P1, P2])
seq3 = generate(g3, 10_000; rng)

# Property-based: intermittent-map latent driver + quantization
g3b = IntermittentMapSymbols(1.6, [:a, :b, :c], [0.2, 0.3, 0.5])
seq3b = generate(g3b, 10_000; rng)

# Model-based: Linear-Additive Markov Process
g4 = LAMP(0.5, [:a, :b, :c], [0.2, 0.3, 0.5]; d = 500, epsilon = 0.02)
seq4 = generate(g4, 10_000; rng)

# Model-based: scalable dyadic-bucket LAMP approximation
g4b = DyadicLAMP(0.5, [:a, :b, :c], [0.2, 0.3, 0.5]; d = 100_000)
seq4b = generate(g4b, 10_000; rng)

# Model-based: centered additive Markov memory function
g4c = CalibratedAdditiveMarkov(0.5, [:a, :b, :c], [0.2, 0.3, 0.5]; d = 500)
seq4c = generate(g4c, 10_000; rng)

# Model-based: heavy-tailed regime-switching Markov chain
Q = [0.2 0.8; 0.8 0.2]
g5 = OnOffMarkov(1.5, [:a, :b], [P1, P2], Q)
seq5 = generate(g5, 10_000; rng)

# Model-based: Fractal Symbol Sequence  (H = (3-α)/2 = 0.75)
g6 = FSS(1.5, [:a, :b, :c]; rates = [2.0, 3.0, 5.0])
seq6 = generate(g6, 10_000; rng)

# Model-based: Hawkes-style self-exciting symbolic process
g7 = HawkesSymbol(0.6, [:a, :b, :c]; d = 500)
seq7 = generate(g7, 10_000; rng)

# Model-based: copy-mutate symbolic growth
g8 = DuplicationMutation(1.5, ['A', 'C', 'G', 'T']; mutation_probability = 0.02)
seq8 = generate(g8, 10_000; rng)

empirical_marginal(seq1, g1.alphabet)

# Save to INC format with full provenance metadata
save_sequence("seq_pb1.inc", seq1, g1)
```

Use `make_generator(id, alphabet; kwargs...)` for standard cases. The `id` may
be a method identifier such as `:PB1`, a string such as `"MB1c"`, or an exported
type-name alias such as `:SpectralFGN`. Use `method_ids()` to list methods and
`method_info(id)` to inspect defaults, standard cases, and a short description.
The method-specific constructors remain the precise API when you need full
control over transition matrices, excitation matrices, or other scientific
settings.

## Standard Factory API

`make_generator` gives all implemented methods one lightweight construction
path for common examples and smoke tests:

```julia
g = make_generator(:MB1c, [:a, :b]; beta = 0.5, d = 500)
seq = generate(g, 1000; rng = StableRNG(7))
```

The factory accepts the following method IDs:

| ID | Type | Standard cases |
|----|------|----------------|
| `:PB1` | `SpectralFGN` | `:standard` |
| `:PB2` | `LGCM` | `:standard` |
| `:PB3` | `WaveletMarkov` | `:persistent_regimes`, `:iid_regimes` |
| `:PB4` | `IntermittentMapSymbols` | `:standard` |
| `:MB1a` | `LAMP` | `:repeat`, `:iid` |
| `:MB1b` | `DyadicLAMP` | `:repeat`, `:iid` |
| `:MB1c` | `CalibratedAdditiveMarkov` | `:standard`, `:iid` |
| `:MB2` | `OnOffMarkov` | `:persistent_regimes`, `:iid_regimes` |
| `:MB3` | `FSS` | `:standard` |
| `:MB4` | `HawkesSymbol` | `:identity_excitation` |
| `:MB5` | `DuplicationMutation` | `:standard` |

Factory defaults are intentionally conservative standard cases, not scientific
calibrations. For example, the PB3 and MB2 persistent-regime defaults create one
iid and one repeat-biased Markov regime so symbol-level diagnostics can see the
latent regime process. Inspect the exact defaults with:

```julia
method_info(:MB2).defaults
```

## Generators

| ID  | Type          | LRD mechanism                          | Short-range control        | Complexity      |
|-----|---------------|----------------------------------------|----------------------------|-----------------|
| PB1 | `SpectralFGN` | Spectral $1/f^\alpha$ shaping          | Poor (set by quantization) | $O(n \log n)$  |
| PB2 | `LGCM`        | Latent fGn streams + argmax            | Offset-calibrated marginals | $O(n \cdot k \cdot I)$ |
| PB3 | `WaveletMarkov` | Spectral or Haar latent driver + Markov regimes | Per-regime Markov matrices | $O(n \log n + n \cdot k)$ |
| PB4 | `IntermittentMapSymbols` | Intermittent-map latent driver | Poor (set by quantization) | $O(n \log n)$ |
| MB1a | `LAMP`       | Exact power-law history weights        | History-weighted transition matrix | $O(n \cdot \min(d,n))$ |
| MB1b | `DyadicLAMP` | Dyadic approximation to power-law history | History-weighted transition matrix | $O(n k \log n \log \min(d,n))$ |
| MB1c | `CalibratedAdditiveMarkov` | Centered additive memory function | Symbol recurrence through additive memory | $O(n \cdot \min(d,n))$ |
| MB2 | `OnOffMarkov` | Heavy-tailed regime sojourns           | Per-regime Markov matrices | $O(n \cdot k)$ |
| MB3 | `FSS`         | Pareto renewal process per symbol      | Poor (independent streams) | $O(n \cdot k)$ |
| MB4 | `HawkesSymbol` | Power-law self/cross-excitation       | Excitation matrix | $O(n \cdot k \cdot \min(d,n))$ |
| MB5 | `DuplicationMutation` | Power-law lag copy/mutate growth | Growth-induced copy structure | $O(n \log d + d)$ |

### LRD parameters

Each generator exposes either a direct Hurst parameter $H \in (1/2, 1)$ or a
tail/decay parameter with a nominal relationship to $H$:

| Type          | Input parameter  | Relationship to $H$    |
|---------------|------------------|------------------------|
| `SpectralFGN` | `H`              | direct                 |
| `LGCM`        | `H`              | direct                 |
| `WaveletMarkov` | `H`            | latent-driver target   |
| `IntermittentMapSymbols` | `z` | latent intermittency strength |
| `LAMP`        | `beta` ($\beta$) | $H = (2 - \beta) / 2$ |
| `DyadicLAMP`  | `beta` ($\beta$) | finite dyadic approximation to $H = (2 - \beta) / 2$ |
| `CalibratedAdditiveMarkov` | `beta` ($\beta$) | finite additive memory-function decay |
| `OnOffMarkov` | `alpha` ($\alpha$) | nominal $H = (3 - \alpha) / 2$ |
| `FSS`         | `alpha` ($\alpha$) | $H = (3 - \alpha) / 2$ |
| `HawkesSymbol` | `beta` ($\beta$) | finite power-law excitation kernel |
| `DuplicationMutation` | `alpha` ($\alpha$) | copy-distance exponent; validation uses `alpha - 1` as an empirical reference slope |

## Controllability

All implemented generators accept an explicit ordered alphabet. Duplicate alphabet
entries are rejected because they make empirical frequency tables ambiguous.

Use `control_capabilities(g)` to inspect the strength and scope of the controls a
generator claims. It distinguishes direct finite-sample or empirical controls from
implied, asymptotic, induced, latent, and nominal behavior.

| Type | Marginal control | Local structure control |
|------|------------------|-------------------------|
| `SpectralFGN` | direct `marginal`; rank binning gives integer counts as close as possible to target | none |
| `LGCM` | direct `marginal`; calibrated latent offsets | none |
| `WaveletMarkov` | aggregate stationary marginal implied by regimes | per-regime bigram matrices |
| `IntermittentMapSymbols` | direct `marginal`; rank binning gives integer counts as close as possible to target | none |
| `LAMP` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control | exact history-weighted transition matrix |
| `DyadicLAMP` | direct `marginal` mixed through `epsilon`; larger `epsilon` improves marginal control | dyadic-bucket approximation to history-weighted transition matrix |
| `CalibratedAdditiveMarkov` | centered additive memory around `marginal`; `strength = 0` is iid | additive memory induces recurrence |
| `OnOffMarkov` | aggregate stationary marginal implied by regimes | per-regime bigram matrices |
| `FSS` | asymptotic `rates / sum(rates)` | none |
| `HawkesSymbol` | baseline distribution reported, output marginal implied by excitation | excitation matrix |
| `DuplicationMutation` | seed/mutation replacement distribution; output marginal shaped by copied history | copy/mutate growth-induced |

First-order local structure can be represented by a validated `MarkovSpec`:

```julia
spec = MarkovSpec([:a, :b], [0.9 0.1; 0.2 0.8])
g1 = WaveletMarkov(0.8, [spec, spec])
g2 = OnOffMarkov(1.5, [spec, spec], [0.2 0.8; 0.8 0.2])

control_capabilities(g1).bigram
```

`WaveletMarkov` and `OnOffMarkov` provide `:per_regime` bigram control. Their
aggregate observed bigrams depend on regime differences and switching behavior. If
all regimes use the same `MarkovSpec`, that common transition matrix is also the
unambiguous aggregate target.

`MarkovSpec` is the current concrete `LocalStructureSpec`, and
`local_structure_order(spec)` returns its order. This is the intended extension
point for future sparse higher-order controls. S5.jl currently provides
`empirical_trigram` for diagnostics, but it does not expose a trigram-control
specification.

For `WaveletMarkov` and `OnOffMarkov`, one-hot symbol diagnostics need regimes
with different stationary symbol distributions. If each regime has the same
stationary marginal, the latent regime process can be long-memory while the
symbol-level ACF and spectrum look nearly short-memory.
`WaveletMarkov` defaults to `driver = :spectral`, which rank-bins an approximate
spectral fGn latent series into regimes. The legacy `driver = :haar` cascade is
retained for comparison in validation studies.

Reproducible controllability studies live in `validation/`, for example:

```julia
julia --project=. validation/marginal_control.jl
julia --project=. validation/local_structure.jl
julia --project=. validation/lrd_method_diagnostics.jl
julia --project=validation validation/longmemory_comparison.jl
```

The LRD diagnostic script writes generated sequences and summary tables as INC
files under `validation/results/lrd_diagnostics/`, and writes log-log SVG plots of
log-binned one-hot autocorrelation and power-spectrum summaries.
The symbolic-to-numeric transformation is formalized in
`validation/lrd_symbol_diagnostics.jl`: each symbol is converted to a centered
one-hot indicator series before autocorrelation, autocovariance, or periodogram
calculations. `validation/longmemory_comparison.jl` compares those helpers with
LongMemory.jl's `autocovariance`, `autocorrelation`, and `periodogram`, including
the documented lag-zero and angular-frequency adaptations.
Autocorrelation SVG plots include dashed vertical interpretation limits: a
finite-sample `n / 10` lag limit, and explicit generator limits where they exist,
such as `LAMP.d`. Power-spectrum plots show the same scales as reciprocal
frequencies. Autocorrelation plots also include a gray dashed nominal
power-law reference line with slope `lag^(-beta)`, anchored to the first positive
plotted autocorrelation value. Power-spectrum plots include the corresponding
gray dashed low-frequency reference with slope `frequency^(beta - 1)`.

See `VALIDATION_POLICY.md` for the validation tiers. Fast tests run through the
package test suite. Longer validation studies and larger benchmark runs are manual
or flag-controlled.

## Benchmarking

Benchmarking uses a separate environment under `benchmark/`:

```julia
julia --project=benchmark benchmark/benchmarks.jl
S5_BENCHMARK_LARGE=true julia --project=benchmark benchmark/benchmarks.jl
```

The default suite covers all implemented generators across moderate sequence
lengths and alphabet sizes. The large suite adds longer runs and should be treated
as machine-specific performance evidence rather than a correctness test.

## Motivations

S5.jl is intended for several related research uses:

- testing NN-LRD estimators on sequences with known generator settings;
- probing excess entropy rate, entropy-rate convergence, and other
  information-theoretic summaries of long-memory symbolic data;
- creating non-language but language-like long-context sequences for training and
  stress-testing LLM-style neural networks;
- separating local bigram/trigram competence from genuine long-context modelling;
- studying anomaly detection, change detection, and synthetic privacy-preserving
  categorical data with realistic burstiness.

## INC output format

All generated sequences should be saved with `save_sequence`, which writes the
data and full provenance metadata (package version, generator type and
parameters, creation date) in the
[INC format](https://github.com/mroughan/IncCSV.jl).

```julia
save_sequence("output.inc", seq, g)
```

The resulting file has this structure:

```
---
title = S5.jl synthetic LRD symbol sequence
package = S5
version = "0.1.0"
generator = SpectralFGN
method = PB1
n = 10000
created = "2026-06-04"
[generator_params]
H = "0.8"
alphabet_size = 3
alphabet = "a,b,c"
marginal = "0.33333333,0.33333333,0.33333334"
[columns]
index = time index (1-based)
symbol = generated symbol from the alphabet
---
index,symbol
1,b
2,a
…
```

## Background

Long-Range Dependence (LRD) means that the large-scale statistical structure
of a sequence is as important as its short-range structure. Formally, a
sequence is LRD if its autocovariance function (ACF) decays as a power law:

$$\gamma_k \sim c_\gamma |k|^{-\beta}, \quad \beta \in (0, 1),$$

with Hurst parameter $H = (2-\beta)/2 \in (1/2, 1)$.

LRD is ubiquitous in human-generated data (natural language, Internet traffic,
genomics, social media), yet almost all LRD synthesis tools target numerical
data. S5.jl fills that gap for **symbol sequences** — data taking values in a
finite, unordered alphabet.

## References

- Paxson, V. (1997). Fast, approximate synthesis of fractional Gaussian noise. *CCR* 27.
- Li, W. (1991). Expansion-modification systems: a model for spatial 1/f spectra. *Physical Review A* 43.
- Li, W., Marr, T. G., & Kaneko, K. (1994). Understanding long-range correlations in DNA sequences. *Physica D* 75.
- Melnyk, S. S., Usatenko, O. V., & Yampol'skii, V. A. (2006). Memory functions of the additive Markov chains. *Physica A* 361.
- Mayzelis, Z. A., Apostolov, S. S., Melnyk, S. S., Usatenko, O. V., & Yampol'skii, V. A. (2006). Additive N-step Markov chains as prototype model of symbolic stochastic dynamical systems with long-range correlations.
- Kumar, R., Raghu, M., Sarlós, T., & Tomkins, A. (2017). Linear additive Markov processes. *WWW '17*.
- Lowen, S. B. & Teich, M. C. (1995). Fractal stochastic point processes. *Fractals* 3(1).
- Provata, A., & Beck, C. (2012). Coupled intermittent maps modelling the statistics of genomic sequences: a network approach. arXiv:1205.2249.
- Pipiras, V. & Taqqu, M. S. (2017). *Long-Range Dependence and Self-Similarity*. Cambridge UP.
