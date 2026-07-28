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
- **Type-stable:** instruction sequences are tuples of concrete types; results
  are `NamedTuple`s that the Julia compiler can specialize on.
- **Modular backend:** the core DSL is backend-agnostic; the current physics
  backend is `FourVectors.jl`.
- **Unified execution:** [`apply_decay_instruction`](@ref) is the single public
  entry point for one instruction, a tuple of instructions, or a
  [`CompositeInstruction`](@ref).

## Getting started

```julia
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

```julia
composite = CompositeInstruction((
    ToHelicityFrame((1, 2, 3)),
    PlaneAlign(4, 5),
))

final_objs, results = apply_decay_instruction(composite, objs)
```

For path comparisons, use [`compare_instruction_paths`](@ref) and extract
relative Wigner angles with [`wigner_zyz`](@ref).

## Index specification

Instructions that take indices accept:

- a single `Int` (one four-vector),
- a `Tuple` or `Vector{Int}` (sum of the selected four-vectors), or
- varargs such as `ToHelicityFrame(1, 2)` on supported constructors.

Negative indices use the corresponding four-vector with a minus sign. This is
useful for axis alignment and for comparing equivalent path conventions:

```julia
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

### Carried helicity axes

[`ToHelicityFrame`](@ref) realigns the frame at every step, so the helicity ẑ is
always `(0,0,1)` and nothing has to be remembered. The alternative — used by
TF-PWA — is to boost *without* realigning and carry the axes along as data:

```julia
ToHelicityFrame  ==  ToRestFrame |> Rz(-ϕ) |> Ry(-θ)
```

Those two rotations are the entire difference. [`ToRestFrame`](@ref) is the pure
boost; the alignment rotation it declines to apply is instead recorded in a
[`HelicityAxes`](@ref) marker, which occupies a slot of `objs` added by
[`with_helicity_axes`](@ref) and is inert under pure boosts. The same program
written the other way:

```julia
program = (
    PlantLabAxes(5),                # ẑ=(0,0,1), x̂=(1,0,0)   [axes planted]
    ToRestFrame((1, 2, 3, 4)),      # pure boost              [total rest frame]
    TransportAxes(5, (1, 2, 3)),    # carry axes onto (1,2,3) [one vertex]
    ToRestFrame((1, 2, 3)),         # pure boost              [(1,2,3) rest frame]
    MeasureEulerZXZ(:v, (1, 2), 5), # (α, β) against carried axes
)
final_objs, results = apply_decay_instruction(program, with_helicity_axes(objs))
```

- [`ToRestFrame`](@ref): pure boost, no realignment — see [`boost_to_rest`](@ref).
- [`PlantLabAxes`](@ref): write the starting ẑ, x̂ convention into a slot.
- [`TransportAxes`](@ref): carry `(ẑ, x̂)` across one decay vertex.
- [`MeasureEulerZXZ`](@ref): [`MeasureCosThetaPhi`](@ref) against the carried
  axes; `α` plays the role of `ϕ` and `β` of `θ`, via [`euler_zxz`](@ref).

Both routes give the same angles. The carried-axis form is useful when comparing
against frameworks that build angles this way, and when the parent is at rest —
a pure boost with `γ = 1` is exactly the identity, whereas realigning applies the
undefined `ϕ`, `θ` of a null three-momentum.

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
