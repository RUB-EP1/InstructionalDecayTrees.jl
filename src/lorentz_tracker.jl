"""
    LorentzTracker([T=Float64])

Tracks an accumulated Lorentz transformation matrix `Λ` acting on column vectors
ordered as `(px, py, pz, E)`.
Also carries the corresponding 2x2 SU(2) matrix `U` for phase-aware tracking.
"""
struct LorentzTracker{T<:Real,TM4<:AbstractMatrix{T},TM2<:AbstractMatrix{Complex{T}}}
    Λ::TM4
    U::TM2
end

function LorentzTracker(::Type{T}=Float64) where {T<:Real}
    return LorentzTracker{T,Matrix{T},Matrix{Complex{T}}}(
        Matrix{T}(I, 4, 4),
        Matrix{Complex{T}}(I, 2, 2),
    )
end

Base.inv(t::LorentzTracker) = LorentzTracker(inv(t.Λ), inv(t.U))
Base.:*(a::LorentzTracker, b::LorentzTracker) = LorentzTracker(a.Λ * b.Λ, a.U * b.U)

"""
    relative_tracker(reference, other)

Relative transform that maps vectors expressed in the reference path frame
to the corresponding vectors in the other path frame:

`Δ = other * inv(reference)`
"""
relative_tracker(reference::LorentzTracker, other::LorentzTracker) = other * inv(reference)

"""
    decode_lorentz_helicity(t; atol=10*eps(T))

Decode full Lorentz parameters in helicity convention from a tracked transform.
Returns `(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf)`.

The parameters are decoded from the 4×4 matrix `Λ` alone, so the result is
branch-blind: trackers with `U` and `-U` decode identically (`ψ_rf` is
determined modulo 2π, not 4π). For spinor-branch-exact Wigner angles of a
pure-rotation tracker, use [`wigner_zyz`](@ref).

Default `atol` is `10 * eps(T)` where `T = eltype(t.Λ)`.
"""
function decode_lorentz_helicity(t::LorentzTracker; atol::Real=_default_atol(t.Λ))
    d = _decode_lorentz_helicity_zyz_xyze(t.Λ; atol = atol)
    return (
        ϕ = d.ϕ,
        θ = d.θ,
        ξ = d.ξ,
        ϕ_rf = d.ϕ_rf,
        θ_rf = d.θ_rf,
        ψ_rf = normalize_psi(d.ψ_rf),
    )
end

"""
    wigner_zyz(t; atol=10*eps(T))

Extract active helicity-convention ZYZ Wigner angles `(ϕ, θ, ψ)` of a
pure-rotation tracker, decoded from the tracked SU(2) matrix `t.U`.

Wigner angles are defined for rotations. The relative transform between two
helicity paths to the same particle is a pure rotation by construction; a
tracker that still carries a boost (`ξ` above tolerance) raises an error —
use [`decode_lorentz_helicity`](@ref) for general transforms.

The spinor ±1 branch is intrinsic to the SU(2) decode: `ψ` is returned on
`[-π, 3π)` and the angles rebuild `t.U` exactly (never `-t.U`), so Wigner
`D^j` matrices built from them are phase-exact for half-integer `j`.

Default `atol` is `10 * eps(T)` where `T = eltype(t.Λ)`.

For an instructive comparison of SO(3) (`Λ`) vs SU(2) (`U`) decoders, see
`docs/wigner_su2_so3.qmd`.
"""
function wigner_zyz(t::LorentzTracker; atol::Real=_default_atol(t.Λ))
    b = _decode_boost_xyze(t.Λ; atol = atol)
    abs(b.ξ) < atol || throw(ArgumentError(
        "wigner_zyz decodes Wigner angles of a pure rotation, but this tracker " *
        "carries a boost (ξ ≈ $(b.ξ)). Use decode_lorentz_helicity for general transforms."))
    return _decode_rotation_zyz_su2(t.U; atol = atol)
end

"""
    _wigner_zyz_so3(t; atol=10*eps(T))

Internal: ZYZ angles from the spatial SO(3) block of `Λ` only (no SU(2) branch).
Used by the Wigner-angle tutorial; prefer [`wigner_zyz`](@ref).
"""
function _wigner_zyz_so3(t::LorentzTracker; atol::Real=_default_atol(t.Λ))
    d = _decode_lorentz_helicity_zyz_xyze(t.Λ; atol = atol)
    return (ϕ = d.ϕ_rf, θ = d.θ_rf, ψ = normalize_psi(d.ψ_rf))
end

"""
    _wigner_zyz_flip(t; atol=10*eps(T))

Legacy implementation of [`wigner_zyz`](@ref), kept only as an independent
reference for tests — do not use in new code.

Decodes the angles from the SO(3) block of `Λ` (branch-blind), rebuilds the
SU(2) matrix from them, and if the rebuild matches `-t.U` repairs the spinor
branch by shifting `ψ` by 2π. `wigner_zyz` reads the branch directly off `t.U`
and never needs the repair step. On pure rotations the two return the same
triple, except exactly at the `ϕ = ±π` seam, where they may pick different —
equally faithful — representatives (`(ϕ, ψ)` vs `(ϕ ∓ 2π, ψ ± 2π)`).
"""
function _wigner_zyz_flip(t::LorentzTracker; atol::Real=_default_atol(t.Λ))
    d = _decode_lorentz_helicity_zyz_xyze(t.Λ; atol = atol)
    ψ_rf = d.ψ_rf
    if abs(d.ξ) < atol
        U_pred = _build_su2(d.ϕ, d.θ, d.ξ, d.ϕ_rf, d.θ_rf, d.ψ_rf)
        err_plus = sum(abs2, U_pred .- t.U)
        err_minus = sum(abs2, U_pred .+ t.U)
        ψ_rf = err_minus + atol < err_plus ? d.ψ_rf + 2π : d.ψ_rf
    end
    return (ϕ = d.ϕ_rf, θ = d.θ_rf, ψ = normalize_psi(ψ_rf))
end
