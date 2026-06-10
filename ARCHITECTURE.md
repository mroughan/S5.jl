# S5.jl Architecture

## Purpose

S5.jl synthesizes long-range-dependent (LRD) sequences over finite, ordered
alphabets. It is the synthesis component of the research program described in
the project README and grant material. The package provides reproducible
generators with explicit control contracts and provenance, rather than claiming
that every generator controls every statistical property.

This document defines the stable development pathway for S5.jl. The current
research backlog remains in `TODO.md`.

## Goals

- Provide a small, consistent Julia API for generating symbolic LRD sequences.
- Preserve six clearly identified synthesis methods: PB1-PB3 and MB1-MB3.
- State and test what each generator controls, including its alphabet, marginal
  distribution, local structure, and nominal LRD parameter.
- Make stochastic results reproducible when callers supply an explicit random
  number generator.
- Preserve generated-data provenance through `save_sequence` and the INC format.
- Support research use without making unverified scientific claims.
- Maintain good single-thread performance and document algorithmic complexity.
- Keep the package usable on the minimum supported Julia version in
  `Project.toml`.

## Non-Goals

- S5.jl is not an LRD estimator package. Spectral, wavelet, Whittle,
  recurrence-time, Hill, and count-variance estimators belong elsewhere.
- S5.jl does not promise exact recovery of a requested Hurst parameter from a
  finite sequence or from a particular estimator.
- S5.jl does not provide one universal control interface that suggests all
  generators can prescribe arbitrary marginals, bigrams, or trigrams.
- S5.jl is not a general time-series, language-modeling, plotting, or data-frame
  toolkit.
- The fast test suite is not a substitute for Monte Carlo validation or
  scientific review.
- Generated validation sequences are not source artifacts and should not be
  committed unless there is a documented reason to preserve them.

## Philosophy

### Honest Control Contracts

Each generator must document the properties it directly controls, the
properties it only targets asymptotically or empirically, and the properties it
cannot control. Approximate methods must remain labeled approximate. Nominal
relationships between model parameters and LRD measures must not be presented
as finite-sample guarantees.

### Reproducibility Is Part of the API

All stochastic generation must accept `rng::AbstractRNG`. Tests and validation
studies should use `StableRNGs.StableRNG` where reproducibility across Julia
sessions matters. Saved sequences must include sufficient metadata to identify
their generator and parameters.

### Research Evidence Is Separate From Package Correctness

Unit tests establish deterministic contracts, input validation, reproducibility,
and bounded statistical sanity checks. Reproducible studies under `validation/`
provide broader empirical evidence. New scientific claims require a cited
source or a clearly labeled project hypothesis, plus an appropriate validation
study.

### Prefer Small, Explicit Interfaces

The common interface should remain small. Add shared abstractions only when
multiple generators have the same meaningful contract. Do not erase important
method differences merely to make APIs look uniform.

## System Boundaries

### Public Package Layer

- `src/S5.jl` defines the module, exports, dependency imports, and include order.
- `src/interface.jl` defines `LRDGenerator` and the common `generate` interface.
- `src/controls.jl` reports each generator's declared target marginal.
- `src/io.jl` owns stable sequence serialization and provenance metadata.
- `src/utils.jl` contains shared validation, sampling, quantization, and
  empirical-control helpers.

The exported names in `src/S5.jl`, constructor signatures, generated sequence
element types, accepted parameter domains, and INC metadata fields are public
contracts. Changes to them require explicit compatibility consideration,
documentation, tests, and a changelog entry.

### Generator Layer

Each synthesis method lives in one file and is responsible for its own
constructor validation and generation algorithm:

| Family | ID | Type | File |
|---|---|---|---|
| Property-based | PB1 | `SpectralFGN` | `src/pb1.jl` |
| Property-based | PB2 | `LGCM` | `src/pb2.jl` |
| Property-based | PB3 | `WaveletMarkov` | `src/pb3.jl` |
| Model-based | MB1 | `LAMP` | `src/mb1.jl` |
| Model-based | MB2 | `OnOffMarkov` | `src/mb2.jl` |
| Model-based | MB3 | `FSS` | `src/mb3.jl` |

Generator implementations may use private shared helpers, but a generator's
scientific mechanism should remain legible in its own file. A new generator
must have a documented mechanism, control contract, complexity statement,
references, focused tests, validation evidence, and provenance metadata.

### Evidence And Communication Layers

- `test/` contains fast automated contract and regression tests.
- `validation/` contains reproducible, potentially expensive simulation studies.
- `docs/` contains the Documenter.jl user guide and API reference.
- `README.md` is the concise project overview and must retain the AI disclosure.
- `paper/` and `background/` contain research context; package behavior must not
  silently drift to match prose that has not been implemented and tested.
- `TODO.md` records the research and engineering backlog. It is not a guarantee.

## Stable Contracts

### Generator Contract

Every concrete `LRDGenerator` must:

1. Validate parameters at construction time and reject invalid probability
   vectors, alphabets, matrices, and parameter ranges with `ArgumentError`.
2. Implement `generate(g, n; rng = Random.default_rng())`.
3. Return exactly `n` values from the supplied alphabet, preserving its element
   type.
4. Produce repeatable output when invoked with equivalent explicit RNGs.
5. Document its mechanism, control limitations, parameter domains, complexity,
   and scientific references.
6. Implement `target_marginal(g)` when a meaningful declared target exists.
7. Support provenance metadata in `save_sequence`.

### Compatibility Contract

- Follow semantic versioning when changing public behavior.
- Prefer additive changes during the `0.x` research phase.
- Treat exported names, documented constructors, parameter meanings, and output
  metadata as compatibility-sensitive.
- Do not change seeded output casually. When an algorithmic correction requires
  it, document the reason and update reproducibility tests deliberately.
- Avoid new runtime dependencies unless they materially improve correctness,
  reproducibility, or performance.

## Development Pathway

Use this sequence for changes:

1. Define the user-facing or scientific claim and identify whether it belongs in
   the core package, validation studies, or a separate estimator package.
2. Cite the source for scientific or important engineering decisions. Clearly
   label uncertainty and project-specific hypotheses.
3. Implement the narrowest change consistent with existing interfaces.
4. Add or update focused tests for each changed contract.
5. Add or update a validation study when the claim is statistical and cannot be
   established reliably by a fast unit test.
6. Update docstrings, user documentation, `TODO.md` when relevant, and
   `CHANGELOG.md`.
7. Run the package tests on the active Julia environment. Build the
   Documenter.jl docs when exported APIs or docs change.
8. Review performance and allocations for changes on generation hot paths.

## Testing And Validation

- Keep `test/runtests.jl` deterministic and suitable for continuous integration.
- Prefer exact assertions for API contracts and validation errors.
- Use bounded statistical tests sparingly; fix seeds and allow defensible
  tolerances.
- Put broad grids, many replicates, plots, and exploratory diagnostics in
  `validation/`.
- Keep raw or generated data untouched when it is retained. Record provenance
  for external data and reference material.
- Treat validation failures as evidence to investigate, not merely thresholds to
  loosen.

## Coding Guidelines

- Follow the Julia style guide and performance tips linked below.
- Use US English in code and documentation.
- Include docstrings with testable examples for exported types and methods.
- Prefer type-generic public behavior while avoiding avoidable allocations in
  generation loops.
- Keep shared validation and repeated mechanics in focused helpers.
- Preserve the distinction between public exported APIs and private helpers.
- Use mathematical notation in Markdown rather than Unicode approximations when
  explaining formulas.

## Decision Priorities

When goals conflict, use this order:

1. Scientific correctness and honest uncertainty.
2. Reproducibility and provenance.
3. Public API stability.
4. Correctness and input validation.
5. Clear control contracts and documentation.
6. Single-thread performance and memory use.
7. Convenience and feature breadth.

## References

- Project scope and scientific context: `README.md`, `TODO.md`,
  `paper/main.tex`, and `background/Discovery_Grant_application_in_2023__minus_Lewis_.pdf`.
- Julia contributors, [Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/).
- Julia contributors, [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/).
- Semantic Versioning, [Semantic Versioning 2.0.0](https://semver.org/).
- Keep a Changelog, [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/).
- Local, provenance-labeled snapshots of these external sources:
  `references/README.md`.
