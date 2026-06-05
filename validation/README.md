# S5.jl Validation Studies

This folder contains reproducible simulation studies for generator controllability.
These are not LRD estimators. They test whether a generator respects user-facing
controls such as alphabet membership and target marginal probabilities.

The scripts use `StableRNGs.StableRNG` so results are reproducible across Julia
sessions and package updates.

## Marginal Control

Run from the package root:

```julia
julia --project=. validation/marginal_control.jl
```

The script prints aggregate total-variation and maximum absolute marginal errors for
`SpectralFGN`, `LGCM`, `WaveletMarkov`, `LAMP`, `OnOffMarkov`, and `FSS` across
a small grid of sequence lengths, alphabet sizes, and marginal distributions.

`LGCM` is more expensive than the other methods because it calibrates latent
offsets over an `n × k` matrix. Increase `replicates`, `ns`, or
`calibration_iters` manually when running longer studies.

Keep generated data out of the repository unless a result table is intentionally being
tracked.

## LRD Method Diagnostics

Run from the package root:

```julia
julia --project=. validation/lrd_method_diagnostics.jl
```

This generates 30 sequences of length 100,000 for each implemented method using the
alphabet `{A,B,C,D,E}`. Sequences are saved as INC files under
`validation/results/lrd_diagnostics/sequences/`.

For PB3 (`WaveletMarkov`) and MB2 (`OnOffMarkov`), the diagnostic uses five
observable regimes whose stationary distributions are each biased toward one
alphabet symbol. The regimes are balanced overall, but the symbol-level
one-hot ACF and spectrum can see the long-memory regime process. If all regimes
share the same stationary marginal, these diagnostics can look short-memory even
when the latent regime process has long-range structure.

The script computes the one-hot symbol autocorrelation and power spectrum for each
sequence, averages across symbols and replicates, and writes:

- `average_autocorrelation.inc`: full signed averaged autocorrelation by lag;
- `average_power_spectrum.inc`: full averaged periodogram by Fourier frequency;
- `plot_autocorrelation_logbins.inc`: positive log-binned values used for plotting;
- `plot_power_spectrum_logbins.inc`: log-binned spectrum values used for plotting.

Log-log SVG plots are written under `validation/results/lrd_diagnostics/plots/`.
