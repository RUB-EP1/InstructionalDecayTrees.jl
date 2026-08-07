# InstructionalDecayTrees.jl

InstructionalDecayTrees provides declarative decay-angle programs: boost and
rotate four-vectors through instruction paths, measure kinematic variables, and
track accumulated Lorentz transforms for cross-checks between topologies.

## Installation

```julia
using Pkg
Pkg.add("InstructionalDecayTrees")
```

The package depends on [`FourVectors.jl`](https://github.com/mmikhasenko/FourVectors.jl),
which is registered in the General registry and is installed automatically.

## Overview

- **Declarative:** describe boosts, rotations, and angle measurements as an
  instruction path instead of hand-written matrix algebra.
- **Type-stable state flow:** instruction sequences and tree branches are tuples
  of concrete types, allowing Julia to specialize on the program structure.
  Measurements are collected in `NamedTuple`s keyed by the supplied tags.
- **Modular backend:** the core DSL is backend-agnostic; the current physics
  backend is `FourVectors.jl`.
- **Unified execution:** [`apply_decay_instruction`](@ref) is the single public
  entry point for one instruction, a tuple of instructions, or a
  [`CompositeInstruction`](@ref).

## Getting started

```@example getting-started
using InstructionalDecayTrees
using FourVectors

p1 = FourVector(1.0, 0.0, 0.0; M = 0.14)
p2 = FourVector(-1.0, 0.0, 0.0; M = 0.14)
objs = (p1, p2)

program = (
    ToHelicityFrame((1, 2)),
    MeasurePolar(:theta1, 1),
    MeasureInvariant(:m12, (1, 2)),
)

final_objs, results = apply_decay_instruction(program, objs)
```

Tuples are the most convenient way to write a program. You can also wrap a
reusable sequence in a [`CompositeInstruction`](@ref):

```@example getting-started
composite = CompositeInstruction((
    ToHelicityFrame((1, 2)),
    MeasureCosThetaPhi(:angles1, 1),
))

final_objs, results = apply_decay_instruction(composite, objs)
```

For path comparisons, use [`compare_instruction_paths`](@ref) and extract
relative Wigner angles with [`wigner_zyz`](@ref).
Programs with shared prefixes are covered in the
[forked execution tutorial](@ref forked-execution).

## Index specification

Instructions that take indices accept:

- a single `Int` (one four-vector),
- a `Tuple` or `Vector{Int}` (sum of the selected four-vectors), or
- varargs such as `ToHelicityFrame(1, 2)` on supported constructors.

Negative indices use the corresponding four-vector with a minus sign. This is
useful for axis alignment and for comparing equivalent path conventions:

```@example getting-started
ToHelicityFrame(-1)              # uses -objs[1]
PlaneAlign(1, -2)                # z from objs[1], x from -objs[2]
ToHelicityFrame((1, -2, 3))      # objs[1] + (-objs[2]) + objs[3]
```

## Main instructions

### Frame transformations

- [`ToHelicityFrame`](@ref): boost all objects to the rest frame of the sum of
  `indices`.
- [`ToHelicityFrameParticle2`](@ref): same rest frame, but using the
  particle-2 helicity convention (momentum oriented along `-z` before the boost).
- [`PlaneAlign`](@ref): rotate so `z_idx` points along `+z` and `x_idx` lies in
  the `xz` plane.
- [`ToGottfriedJacksonFrame`](@ref): `ToHelicityFrame(system_indices)` followed
  by `PlaneAlign(beam_idx, -target_idx)` for the standard Gottfried–Jackson
  frame.

### Measurements

- [`MeasurePolar`](@ref), [`MeasureSpherical`](@ref), [`MeasureCosThetaPhi`](@ref),
  [`MeasureMassCosThetaPhi`](@ref), and [`MeasureInvariant`](@ref) store
  kinematic quantities in the returned `NamedTuple`.

### Tracked Lorentz execution

For topology cross-checks, run instruction paths on a [`TrackedState`](@ref) to
accumulate a [`LorentzTracker`](@ref):

- [`init_tracked_state`](@ref) starts from identity.
- [`compare_instruction_paths`](@ref) runs reference and alternate paths and
  returns trackers, results, and transformed objects.
- [`decode_lorentz_helicity`](@ref) decodes `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)` in
  helicity convention.
- [`wigner_zyz`](@ref) extracts relative Wigner angles `(ϕ, θ, ψ)`; see the
  [SO(3) vs SU(2) tutorial](@ref wigner) for branch-resolution details.

See also the [API reference](@ref api-reference) for full docstrings.

## Conventions

Rotations in InstructionalDecayTrees are active transformations: they rotate the
four-vectors themselves into the requested frame. Equivalently, a positive
rotation of the objects corresponds to the inverse passive relabeling of the
coordinate axes.

The SO(3) vs SU(2) Wigner-angle walkthrough is generated in CI from the Quarto
source `docs/wigner_su2_so3.qmd` and published as
[Wigner angles: SO(3) vs SU(2)](@ref wigner).

## Developer notes

The repository includes additional material under `docs/`:

- `tracking_and_crosscheck.md` — tracked execution and Python cross-check workflow
- `lorentz_tracker_math.md` — matrix conventions and decode formulas
- `decayangle_crosscheck.md` — reproducibility notes for
  `test/fixtures/decayangle_crosscheck.json`

Cross-check fixtures are validated in `test/crosscheck_json.jl` against
`decayangle` in helicity convention.

To build these docs locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
