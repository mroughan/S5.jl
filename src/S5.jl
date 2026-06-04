"""
    S5

**Self-Similar Symbols Sequence Synthesis.**

S5.jl generates Long-Range Dependent (LRD) sequences of categorical (non-numerical)
symbols for use as ground-truth test data in LRD estimation studies.

# Generators

| ID  | Type          | Mechanism                              |
|-----|---------------|----------------------------------------|
| PB1 | `SpectralFGN` | Spectral fGn synthesis + quantization  |
| MB1 | `LAMP`        | Linear-Additive Markov Process         |
| MB3 | `FSS`         | Fractal Symbol Sequence via FRP/FSNP   |

# Common interface

    generate(g, n; rng = Random.default_rng()) -> Vector
    save_sequence(filepath, seq, g; created)   -> filepath

# References

Roughan, M. & Willinger, W. (2023). Analysis and Synthesis of Long-Range Structure
in Non-Numerical Time Series. ARC Discovery Grant proposal.
"""
module S5

using FFTW: ifft
using Dates: today
using Random
using Statistics: mean, std, quantile
import IncCSV

export LRDGenerator, generate, save_sequence
export SpectralFGN, LAMP, FSS

include("interface.jl")
include("utils.jl")
include("pb1.jl")
include("mb1.jl")
include("mb3.jl")
include("io.jl")

end # module S5
