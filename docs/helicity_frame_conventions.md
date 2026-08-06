# Particle-2 helicity frames for a subchannel resonance


Consider the two orderings

``` text
T1:   ((((4, 5), 3), 2), 1)
T2:   (1, ((3, (4, 5)), 2))
```

Both contain the chain of rest frames

``` text
12345 → 2345 → 345 → 45,
```

but their daughter ordering differs. Here **particle 1 at a vertex**
means the first daughter of that vertex, not necessarily the external
particle with label `1`. A decay-angle program must measure that first
daughter before it descends to the next resonance frame.

## Generate one event

Generate five on-shell four-vectors independently. A fixed seed makes
the tutorial reproducible. Each complete program first applies the same
transformation to the `12345` center-of-mass frame; only subsequent
instructions differ.

``` julia
using InstructionalDecayTrees
using FourVectors
using LinearAlgebra
using Random

rng = MersenneTwister(20260806)
masses = (0.20, 0.30, 0.25, 0.18, 0.22)

lab_objects = Tuple(begin
    momentum = 0.6 .* randn(rng, 3)
    FourVector(momentum[1], momentum[2], momentum[3]; M = mass)
end for mass in masses)

root_to_cmf = ToHelicityFrame((1, 2, 3, 4, 5))
```

## Measure, then descend

The programs spell out the repeated operation explicitly:

1.  arrive in a resonance rest frame;
2.  measure the angles of that vertex’s first daughter;
3.  transform to the next resonance rest frame.

In T1 the successive first daughters are `(2345)`, `(345)`, `(45)`, and
`4`. In T2 they are `1`, `(345)`, `3`, and `4`.

``` julia
T1 = (
    root_to_cmf,
    MeasureSpherical(:theta_12345, :phi_12345, (2, 3, 4, 5)),
    ToHelicityFrame((2, 3, 4, 5)),
    MeasureSpherical(:theta_2345, :phi_2345, (3, 4, 5)),
    ToHelicityFrame((3, 4, 5)),
    MeasureSpherical(:theta_345, :phi_345, (4, 5)),
    ToHelicityFrame((4, 5)),
    MeasureSpherical(:theta_45, :phi_45, 4),
)

T2_particle2 = (
    root_to_cmf,
    MeasureSpherical(:theta_12345, :phi_12345, 1),
    ToHelicityFrameParticle2((2, 3, 4, 5)),
    MeasureSpherical(:theta_2345, :phi_2345, (3, 4, 5)),
    ToHelicityFrame((3, 4, 5)),
    MeasureSpherical(:theta_345, :phi_345, 3),
    ToHelicityFrameParticle2((4, 5)),
    MeasureSpherical(:theta_45, :phi_45, 4),
)

T2_helicity = (
    root_to_cmf,
    MeasureSpherical(:theta_12345, :phi_12345, 1),
    ToHelicityFrameParticle2((2, 3, 4, 5)),
    MeasureSpherical(:theta_2345, :phi_2345, (3, 4, 5)),
    ToHelicityFrame((3, 4, 5)),
    MeasureSpherical(:theta_345, :phi_345, 3),
    ToHelicityFrame((4, 5)),
    MeasureSpherical(:theta_45, :phi_45, 4),
)

_, angles_T1 = apply_decay_instruction(T1, lab_objects)
_, angles_T2_particle2 = apply_decay_instruction(T2_particle2, lab_objects)
_, angles_T2_helicity = apply_decay_instruction(T2_helicity, lab_objects)
```

The two T2 programs differ only in how they approach the `(4,5)` frame.
The standard second-daughter convention uses
[`ToHelicityFrameParticle2`](@ref), while the alternate program uses the
`(4,5)` vector directly in [`ToHelicityFrame`](@ref).

## Angles at every vertex

Each table entry names the first daughter being measured and gives
`(θ, ϕ)` in radians.

``` julia
angle_text(x) = string(round(abs(x) < 5e-13 ? 0.0 : x; digits = 3))
angle_pair(angles, frame) = string(
    "θ = ",
    angle_text(getproperty(angles, Symbol(:theta_, frame))),
    "; ϕ = ",
    angle_text(getproperty(angles, Symbol(:phi_, frame))),
)

rows = (
    (frame = :12345, first = ("2345", "1", "1")),
    (frame = :2345, first = ("345", "345", "345")),
    (frame = :345, first = ("45", "3", "3")),
    (frame = :45, first = ("4", "4", "4")),
)

println("| rest frame | T1 | T2 + ToHelicityFrameParticle2 | T2 + ToHelicityFrame |")
println("|---|---:|---:|---:|")
for row in rows
    values = (
        angle_pair(angles_T1, row.frame),
        angle_pair(angles_T2_particle2, row.frame),
        angle_pair(angles_T2_helicity, row.frame),
    )
    println(
        "| ",
        row.frame,
        " | `",
        row.first[1],
        "`: ",
        values[1],
        " | `",
        row.first[2],
        "`: ",
        values[2],
        " | `",
        row.first[3],
        "`: ",
        values[3],
        " |",
    )
end
```

| rest frame | T1 | T2 + ToHelicityFrameParticle2 | T2 + ToHelicityFrame |
|----|---:|---:|---:|
| 12345 | `2345`: θ = 1.163; ϕ = -0.776 | `1`: θ = 1.979; ϕ = 2.366 | `1`: θ = 1.979; ϕ = 2.366 |
| 2345 | `345`: θ = 2.015; ϕ = 2.953 | `345`: θ = 2.015; ϕ = -0.188 | `345`: θ = 2.015; ϕ = -0.188 |
| 345 | `45`: θ = 1.654; ϕ = 3.015 | `3`: θ = 1.487; ϕ = -0.127 | `3`: θ = 1.487; ϕ = -0.127 |
| 45 | `4`: θ = 1.524; ϕ = -1.199 | `4`: θ = 1.524; ϕ = 1.943 | `4`: θ = 1.524; ϕ = -1.199 |

The table exposes three distinct effects:

1.  In the `12345` frame the two topologies measure opposite daughters:
    T1 measures `(2345)`, while T2 measures `1`. Their directions are
    antipodal.
2.  In the `2345` frame all programs measure `(345)`, but the two
    topology conventions use frames whose azimuths differ by `π`. In the
    `345` frame, T1 measures `(45)` while T2 measures its opposite
    daughter `3`.
3.  In the final `45` frame all programs measure `4`. The standard
    particle-2 approach shifts its azimuth by `π`; T2 with
    `ToHelicityFrame((4,5))` agrees with T1.

## Final-frame comparison

Measurements do not alter the tracked frame transformations, so the
complete programs can also be compared directly:

``` julia
particle2_cmp = compare_instruction_paths(T1, T2_particle2, lab_objects)
helicity_cmp = compare_instruction_paths(T1, T2_helicity, lab_objects)

identity4 = Matrix{Float64}(I, 4, 4)
rotate_z_pi = Diagonal([-1.0, -1.0, 1.0, 1.0])

@assert particle2_cmp.relative.Λ ≈ rotate_z_pi atol = 1e-12
@assert helicity_cmp.relative.Λ ≈ identity4 atol = 1e-12

(
    T2_particle2 = round.(particle2_cmp.relative.Λ; digits = 6),
    T2_helicity = round.(helicity_cmp.relative.Λ; digits = 6),
)
```

    (T2_particle2 = [-1.0 -0.0 0.0 0.0; 0.0 -1.0 -0.0 -0.0; 0.0 0.0 1.0 -0.0; -0.0 0.0 -0.0 1.0], T2_helicity = [1.0 0.0 -0.0 -0.0; -0.0 1.0 -0.0 -0.0; 0.0 0.0 1.0 -0.0; -0.0 -0.0 -0.0 1.0])

Thus, relative to T1, T2 with `ToHelicityFrameParticle2((4,5))` leaves
`Rz(π)`, whereas T2 with `ToHelicityFrame((4,5))` leaves the identity.
This statement concerns the spatial Lorentz frame; its SU(2) lift can
additionally carry the usual central sign associated with a `2π` spinor
rotation.
