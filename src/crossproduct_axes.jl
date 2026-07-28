# Cross-product helicity axes.
#
# `ToHelicityFrame` realigns the frame at every step, so the helicity ẑ is
# always (0,0,1) and nothing has to be remembered.  The alternative — used by
# TF-PWA, and the one the `(α, β)` Euler angles of `angle_zx_z_getx` come from —
# is to boost *without* realigning and carry the axes along as data.
#
# The two are the same convention (see `boost_to_rest`), so this file adds a
# second way to write the same program, not a second convention.

"""
    HelicityAxes(ẑ, x̂)

A pair of orthonormal direction markers recording where the helicity ẑ and x̂ of
the current decay step point, in the coordinates of whatever frame the objects
are currently in.

It occupies one slot of `objs` alongside the four-vectors. It is **not** a
four-vector and must never be summed with one — see [`with_helicity_axes`](@ref)
for how to add slots, and [`ToRestFrame`](@ref) for how it transforms.

`HelicityAxes()` is the "not yet planted" marker (`NaN`), so forgetting
[`PlantLabAxes`](@ref) fails loudly instead of silently using stale axes.
"""
struct HelicityAxes{T<:Real}
    ẑ::NTuple{3,T}
    x̂::NTuple{3,T}
end

HelicityAxes{T}() where {T<:Real} =
    HelicityAxes{T}(ntuple(_ -> T(NaN), 3), ntuple(_ -> T(NaN), 3))
HelicityAxes() = HelicityAxes{Float64}()

function Base.show(io::IO, a::HelicityAxes)
    print(io, "HelicityAxes(ẑ=", a.ẑ, ", x̂=", a.x̂, ")")
end

"""
    with_helicity_axes(objs, naxes = 1)

Append `naxes` unplanted [`HelicityAxes`](@ref) slots to `objs`.

Particle indices `1:length(objs)` are untouched, so a program written against
the bare four-vectors keeps working. This is state construction, not physics:
the convention itself is chosen by [`PlantLabAxes`](@ref) inside the program.

```jldoctest
julia> using InstructionalDecayTrees, FourVectors

julia> objs = (FourVector(1.0, 0.0, 0.0; M = 0.14), FourVector(-1.0, 0.0, 0.0; M = 0.14));

julia> length(with_helicity_axes(objs))
3
```
"""
with_helicity_axes(objs, naxes::Int = 1) =
    (objs..., ntuple(_ -> HelicityAxes(), naxes)...)

"""
    helicity_axes_at(objs, slot)

Read the [`HelicityAxes`](@ref) marker stored in slot `slot` of `objs`.
"""
helicity_axes_at(objs, slot::Int) = objs[slot]::HelicityAxes

_set_slot(objs::Tuple, slot::Int, value) =
    ntuple(k -> k == slot ? value : objs[k], length(objs))

# --- three-vector helpers -------------------------------------------------
# Written on NTuple{3} rather than reaching for a vector library: these six
# lines are part of what the reader is meant to check.

_dot3(a::NTuple{3}, b::NTuple{3}) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_norm3(a::NTuple{3}) = sqrt(_dot3(a, a))
_unit3(a::NTuple{3}) = (n = _norm3(a); (a[1] / n, a[2] / n, a[3] / n))
_cross3(a::NTuple{3}, b::NTuple{3}) = (
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1],
)

"""Degeneracy threshold below which a cross product is treated as collinear."""
const AXIS_EPS = 1.0e-12

"""
    cross_unit(a, b)

Unit normal of `a` and `b`, with a reproducible tie-break for the collinear
case: when `a ∥ b` the normal is undefined, so `b` is nudged along `(1,1,1)`
first. Matches TF-PWA `Vector3.cross_unit`.
"""
function cross_unit(a::NTuple{3}, b::NTuple{3})
    c = _cross3(a, b)
    if _norm3(c) < AXIS_EPS
        c = _cross3(a, (b[1] + 1, b[2] + 1, b[3] + 1))
    end
    return _unit3(c)
end

_direction(p::FourVector) = _unit3((p.px, p.py, p.pz))

"""
    euler_zxz(ẑ₁, x̂₁, ẑ₂) -> (α, β, γ)

Euler angles of the rotation `Rz(α) Ry(β) Rz(γ)` carrying the triad
`(x̂₁, ŷ₁, ẑ₁)` onto a triad whose z-axis is `ẑ₂`. `γ ≡ 0` by construction: the
rotation is taken in the `(ẑ₁, ẑ₂)` plane, which is the helicity convention.
Matches TF-PWA `EulerAngle.angle_zx_z_getx`.

`α` and `β` are just an azimuth and a polar angle. The only difference from
`MeasureCosThetaPhi` is that they are measured against the supplied axes rather
than the coordinate axes — when `(ẑ₁, x̂₁) = ((0,0,1), (1,0,0))` they coincide
with `(ϕ, θ)` exactly.
"""
function euler_zxz(ẑ₁::NTuple{3}, x̂₁::NTuple{3}, ẑ₂::NTuple{3})
    u_ẑ₁ = _unit3(ẑ₁)
    u_ẑ₂ = _unit3(ẑ₂)
    ŷ₁ = cross_unit(u_ẑ₁, x̂₁)          # complete the carried triad
    u_x̂₁ = cross_unit(ŷ₁, u_ẑ₁)        # re-orthogonalise x̂₁ against ẑ₁
    ŷ_r = cross_unit(u_ẑ₁, u_ẑ₂)       # normal of the (ẑ₁, ẑ₂) plane
    x̂_r = cross_unit(ŷ_r, u_ẑ₁)        # in-plane companion of ẑ₁
    α = atan(_dot3(x̂_r, ŷ₁), _dot3(x̂_r, u_x̂₁))   # azimuth of the plane about ẑ₁
    β = atan(_dot3(u_ẑ₂, x̂_r), _dot3(u_ẑ₂, u_ẑ₁)) # polar angle of ẑ₂ from ẑ₁
    return (α = α, β = β, γ = 0.0)
end

"""
    _transported_x(ẑ₁, ẑ₂)

The new x̂ after the helicity frame is carried from `ẑ₁` onto `ẑ₂`: rotate about
the `(ẑ₁, ẑ₂)` plane normal, so `x̂₂ = ŷ_r × ẑ₂`. Independent of the old `x̂₁` —
the azimuth is fixed by the decay plane, not by where x̂ used to point.
"""
_transported_x(ẑ₁::NTuple{3}, ẑ₂::NTuple{3}) = cross_unit(cross_unit(ẑ₁, ẑ₂), ẑ₂)

"""
    _momentum(objs, indices)

`get_fourvector` with a guard: axis markers are not four-vectors and summing one
into a momentum is always a bug, so say so rather than failing in `+`.
"""
function _momentum(objs, indices::Tuple{Vararg{Int}})
    for i in indices
        objs[abs(i)] isa FourVector || throw(ArgumentError(
            "slot $(abs(i)) holds $(typeof(objs[abs(i)])), not a FourVector"))
    end
    return get_fourvector(objs, indices)
end

# --- the pure boost -------------------------------------------------------

"""
    boost_to_rest(p, P)

Pure boost of the four-vector `p` into the rest frame of `P`: rotate `P` onto
+ẑ, boost along ẑ, rotate back.

```julia
p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ) |> Ry(θ) |> Rz(ϕ)      # ϕ, θ, γ of P
```

Drop the two trailing rotations and this is `FourVectors.transform_to_cmf`, i.e.
what [`ToHelicityFrame`](@ref) applies. That is the whole difference between the
two frame conventions:

```julia
ToHelicityFrame  ==  ToRestFrame |> Rz(-ϕ) |> Ry(-θ)
```

`ToHelicityFrame` *spends* the alignment rotation, so its ẑ is always `(0,0,1)`;
`ToRestFrame` keeps it, and the price is carrying [`HelicityAxes`](@ref).

Because the rotations cancel exactly, this is also well-behaved when `P` is
already at rest: `γ == 1` makes `Bz` the identity, whatever `polar_angle` and
`azimuthal_angle` return for a null three-momentum. `transform_to_cmf` has no
such protection — it applies those undefined angles.
"""
function boost_to_rest(p::FourVector, P::FourVector)
    θ = polar_angle(P)
    ϕ = azimuthal_angle(P)
    γ = boost_gamma(P)
    return p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ) |> Ry(θ) |> Rz(ϕ)
end

"""
    boost_to_rest(a::HelicityAxes, P)

A direction marker is **inert** under a pure boost.

A pure boost relates two frames whose spatial axes are parallel by definition,
so a marker recording "the previous helicity ẑ pointed there" keeps the same
three numbers. This is what lets the axes live in `objs` next to the
four-vectors without ever leaving the common frame, and it is exactly the
assumption TF-PWA makes when it reuses `set_z[core]` after `cal_chain_boost`.

Under a *rotation* the marker would rotate like any three-vector; no instruction
in this set rotates the frame, which is why it needs no such method.
"""
boost_to_rest(a::HelicityAxes, ::FourVector) = a

# --- instructions ---------------------------------------------------------

"""
    PlantLabAxes(slot)

Write the starting helicity convention into `slot`: `ẑ = (0,0,1)`, `x̂ = (1,0,0)`.

This is where the lab-axis convention enters the program — the cross-product
counterpart of "which frame do we start in?". It is a state write, not a
transform: the four-vectors do not move.
"""
struct PlantLabAxes <: AbstractInstruction
    slot::Int
end

Base.show(io::IO, instr::PlantLabAxes) = print(io, "PlantLabAxes(", instr.slot, ")")

function apply_decay_instruction(instr::PlantLabAxes, objs)
    planted = HelicityAxes((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    return (_set_slot(objs, instr.slot, planted), _empty_instruction_results)
end

"""
    ToRestFrame(indices)
    ToRestFrame(i, j, ...)

Pure boost of everything in `objs` into the rest frame of the sum selected by
`indices`.

The difference from [`ToHelicityFrame`](@ref) is that the frame is *not*
realigned so the boost direction becomes +ẑ — see [`boost_to_rest`](@ref).
Four-vectors move; [`HelicityAxes`](@ref) markers do not. Both end up in the new
frame.

`indices` may be a single integer, a tuple, a vector, or varargs. Negative
indices use the corresponding four-vector with the opposite sign.
"""
struct ToRestFrame{T<:Tuple} <: AbstractInstruction
    indices::T

    function ToRestFrame(indices)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(indices_norm)
    end
end
ToRestFrame(indices::Int...) = ToRestFrame(indices)

function Base.show(io::IO, instr::ToRestFrame)
    print(io, "ToRestFrame(")
    show_indices(io, instr.indices)
    print(io, ")")
end

function apply_decay_instruction(instr::ToRestFrame, objs)
    P = _momentum(objs, instr.indices)
    return (map(o -> boost_to_rest(o, P), objs), _empty_instruction_results)
end

"""
    TransportAxes(slot, along)
    TransportAxes(slot, i, j, ...)

Carry the helicity axes in `slot` across one decay vertex: the new ẑ is the
direction of `sum(objs[along])` in the current frame, and the new x̂ is fixed by
the decay plane.

Run it **while still in the parent rest frame** — `along` must be the daughter's
momentum measured there, as in TF-PWA's `set_z[out] = data[core][out]`. The
following `ToRestFrame(along)` then leaves the fresh axes untouched.

This is the instruction for a vertex that is only on the *path* to the angle you
want. TF-PWA computes `(α, β)` at such vertices too, but
`HelicityDecay.get_D_matrix_term` reads only `data[outs[0]]["ang"]`, so those
never enter an amplitude. Keeping transport separate from measurement makes that
structural: a vertex you pass through gets `TransportAxes`, a vertex you use gets
[`MeasureEulerZXZ`](@ref).
"""
struct TransportAxes{T<:Tuple} <: AbstractInstruction
    slot::Int
    along::T

    function TransportAxes(slot::Int, along)
        along_norm = normalize_indices(along)
        new{typeof(along_norm)}(slot, along_norm)
    end
end
TransportAxes(slot::Int, along::Int...) = TransportAxes(slot, along)

function Base.show(io::IO, instr::TransportAxes)
    print(io, "TransportAxes(", instr.slot, ", ")
    show_indices(io, instr.along)
    print(io, ")")
end

function apply_decay_instruction(instr::TransportAxes, objs)
    old = helicity_axes_at(objs, instr.slot)
    ẑ₂ = _direction(_momentum(objs, instr.along))
    carried = HelicityAxes(ẑ₂, _transported_x(old.ẑ, ẑ₂))
    return (_set_slot(objs, instr.slot, carried), _empty_instruction_results)
end

"""
    MeasureEulerZXZ(tag, indices, slot; branch = 1)

Measure the Euler angles `(α, β, γ=0)` of `sum(objs[indices])` against the axes
carried in `slot`, and store them under `tag` as `(α, β, γ, cosβ)`.

This is [`MeasureCosThetaPhi`](@ref) with the reference axes supplied
explicitly: `α` plays the role of `ϕ` and `β` of `θ`, and when the carried axes
are the coordinate axes the two agree exactly. See [`euler_zxz`](@ref).

`branch` reproduces TF-PWA's per-daughter `bias`: daughter `n` folds `α` into
`[-nπ, (2-n)π)`. Only `branch = 1`, the first-listed daughter, ever reaches a
Wigner D; the keyword exists to document that the other branch is a bookkeeping
convention rather than different physics.
"""
struct MeasureEulerZXZ{T<:Tuple} <: AbstractMeasureInstruction
    tag::Symbol
    indices::T
    slot::Int
    branch::Int

    function MeasureEulerZXZ(tag::Symbol, indices, slot::Int; branch::Int = 1)
        indices_norm = normalize_indices(indices)
        new{typeof(indices_norm)}(tag, indices_norm, slot, branch)
    end
end

function Base.show(io::IO, instr::MeasureEulerZXZ)
    print(io, "MeasureEulerZXZ(")
    show(io, instr.tag)
    print(io, ", ")
    show_indices(io, instr.indices)
    print(io, ", ", instr.slot)
    instr.branch == 1 || print(io, "; branch=", instr.branch)
    print(io, ")")
end

"""Fold daughter `n`'s azimuth into `[-nπ, (2-n)π)`, i.e. TF-PWA's `bias = -nπ`."""
_branch_fold(α, branch::Int) = (bias = -branch * π; mod(α - bias, 2π) + bias)

function apply_decay_instruction(instr::MeasureEulerZXZ, objs)
    ax = helicity_axes_at(objs, instr.slot)
    ẑ₂ = _direction(_momentum(objs, instr.indices))
    ang = euler_zxz(ax.ẑ, ax.x̂, ẑ₂)
    value = (
        α = _branch_fold(ang.α, instr.branch),
        β = ang.β,
        γ = ang.γ,
        cosβ = cos(ang.β),
    )
    return (objs, NamedTuple{(instr.tag,)}((value,)))
end
