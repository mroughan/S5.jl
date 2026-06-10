# Instructions For AI Agents

## Required Reading

Before changing this repository, read:

1. `../GUARDRAILS.md`
2. `ARCHITECTURE.md`
3. `README.md`
4. The relevant source, tests, docs, validation scripts, and `TODO.md`

Follow `../GUARDRAILS.md` if instructions conflict. Never edit that file. Do not
edit `ARCHITECTURE.md` or this file unless the user explicitly asks.

## Project Direction

S5.jl is a focused Julia package for synthesizing LRD symbolic sequences. Keep
estimation methods outside this repository. Preserve honest distinctions between
exact, empirical, asymptotic, latent, and nominal statistical claims.

The stable pathway is:

- keep the common generator API small;
- preserve each synthesis method's scientific identity;
- make controls and limitations explicit;
- require reproducible RNG use and output provenance;
- use fast tests for package contracts and `validation/` for broader scientific
  evidence;
- protect public APIs and document compatibility-sensitive changes.

## Change Workflow

1. Inspect the repository and git status before editing. Work with existing user
   changes and never discard them.
2. Identify the affected public contract and the architecture layer that owns it.
3. Use references for scientific or important engineering decisions. State
   uncertainty clearly and retain downloaded reference provenance when applicable.
4. Make the smallest coherent change using existing Julia and repository patterns.
5. Add a focused test for every changed idea. Add validation evidence for
   statistical claims that are unsuitable for unit tests.
6. Update docstrings and user documentation when behavior changes.
7. Update `CHANGELOG.md` for every change.
8. Run `julia --project=. -e 'using Pkg; Pkg.test()'`.
9. Run `julia --project=docs docs/make.jl` when exported APIs or docs change.
10. Report uncertainties, tests run, and tests not run.

## Implementation Rules

- Follow the Julia style guide and performance tips.
- Use US English throughout.
- Preserve `generate(g, n; rng)` and explicit `AbstractRNG` support.
- Preserve alphabet element types and validate inputs at construction boundaries.
- Include docstrings with testable examples for exported types and methods.
- Keep generation hot paths suitable for strong single-thread performance.
- Add shared helpers only for meaningful repeated behavior.
- Do not introduce an estimator, broad framework, or dependency without showing
  that it belongs within the goals in `ARCHITECTURE.md`.
- Do not present a finite simulation, latent parameter, or nominal relationship as
  proof of LRD behavior.
- Save retained generated output and external data with clear provenance.
- Keep large reproducible generated sequences out of version control.

## Repository Map

- `src/S5.jl`: module, exports, includes
- `src/interface.jl`: common generator interface
- `src/pb*.jl`, `src/mb*.jl`: synthesis implementations
- `src/controls.jl`: declared target marginals
- `src/utils.jl`: shared helpers
- `src/io.jl`: INC output and provenance
- `test/`: fast contract and regression tests
- `validation/`: reproducible empirical studies
- `docs/`: Documenter.jl documentation
- `paper/`, `background/`: research context

## Git And Safety

- Never commit or push.
- Never make irreversible changes or perform large-scale deletion.
- Never create secrets.
- Do not change files outside this repository.
- Preserve the AI disclosure in `README.md`.
- Warn the user clearly if these instructions or `../GUARDRAILS.md` are in danger
  of falling out of context.

