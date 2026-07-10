"""
    AbstractSpinBasis

Supertype for the spin-quantization conventions used by [`SpinState`](@ref).
Concrete subtypes are [`Helicity`](@ref) and [`Canonical`](@ref); they select the
state-preparation element that maps a rest-frame spin state to the carrier
momentum.
"""
abstract type AbstractSpinBasis end

"""
    Helicity()

Helicity basis. The state is prepared as `R(ϕ,θ)·Bz(ξ)|0,λ⟩`, i.e. a rest-frame
state boosted along `+z` and rotated onto the carrier momentum direction. The
projection index of the coefficients is the helicity λ (spin along the momentum).
"""
struct Helicity <: AbstractSpinBasis end

"""
    Canonical()

Canonical basis. The state is prepared as `R(ϕ,θ)·Bz(ξ)·R(ϕ,θ)⁻¹|0,m⟩`, the pure
(rotationless) boost onto the carrier momentum. The projection index of the
coefficients is the spin along the fixed lab `z` axis.
"""
struct Canonical <: AbstractSpinBasis end

"""
    SpinState(basis, p, twos, coeffs)

Spin state of one particle carried by four-vector `p`, expressed in `basis`
([`Helicity`](@ref) or [`Canonical`](@ref)).

The physical state is `|p, 2s, coeffs⟩`: the momentum geometry `(ϕ, θ, ξ)` is
derived from `p`, while `coeffs ∈ ℂ^(2s+1)` hold the spin content together with the
full spinor phase and the `4π` `±1` branch. `twos = 2s` and `length(coeffs) = 2s+1`.
Coefficient `coeffs[k]` corresponds to projection `m = s - (k-1)`, i.e. index `1`
is `m = +s` and index `end` is `m = -s` (see [`projections`](@ref)).

Construct with [`spin_state`](@ref).
"""
struct SpinState{B<:AbstractSpinBasis,T<:Real,F}
    basis::B
    p::F
    twos::Int
    coeffs::Vector{Complex{T}}
end

"""
    spin_state(basis, p, coeffs)
    spin_state(basis, p, twos, m)

Build a [`SpinState`](@ref) on carrier four-vector `p` in `basis`.

The first form takes an explicit coefficient vector; `twos` is inferred as
`length(coeffs) - 1`. The second form builds the projection eigenstate `|m⟩` for
spin `s = twos/2`, where `m` may be integer or half-integer (`2m` must equal an
integer of the same parity as `twos`).
"""
function spin_state(basis::AbstractSpinBasis, p, coeffs::AbstractVector)
    cc = [complex(float(x)) for x in coeffs]
    T = real(eltype(cc))
    return SpinState{typeof(basis),T,typeof(p)}(basis, p, length(cc) - 1, cc)
end

function spin_state(basis::AbstractSpinBasis, p, twos::Integer, m::Real)
    twom = round(Int, 2m)
    (abs(twom) <= twos && iseven(twos - twom)) ||
        throw(ArgumentError("projection 2m=$twom incompatible with twos=$twos"))
    c = zeros(ComplexF64, twos + 1)
    c[(twos - twom) ÷ 2 + 1] = one(ComplexF64)
    return spin_state(basis, p, c)
end

"""
    projections(twos) -> Tuple

Projection values `m` (as `Rational{Int}`) aligned with the coefficient order of a
spin-`twos/2` [`SpinState`](@ref): `(s, s-1, …, -s)`.
"""
projections(twos::Integer) = ntuple(i -> (twos - 2 * (i - 1)) // 2, twos + 1)

# --- internals ---

# 2x2 SL(2,C) preparation elements. `p` supplies (ϕ, θ, ξ). Helicity uses two
# angles (the third ZYZ angle commutes with Bz and lives in `coeffs`); canonical
# is the rotationless boost R·Bz·R⁻¹.
function _prep_su2(::Helicity, p)
    ϕ = azimuthal_angle(p)
    θ = polar_angle(p)
    ξ = acosh(boost_gamma(p))
    return _su2_rz(ϕ) * _su2_ry(θ) * _su2_bz(ξ)
end

function _prep_su2(::Canonical, p)
    ϕ = azimuthal_angle(p)
    θ = polar_angle(p)
    ξ = acosh(boost_gamma(p))
    return _su2_rz(ϕ) * _su2_ry(θ) * _su2_bz(ξ) * _su2_ry(-θ) * _su2_rz(-ϕ)
end

# (2s+1)-dim representation of a 2x2 matrix `w` via the Schwinger symmetric-power
# construction. Homomorphism for any 2x2 matrix, and `_wignerD(1, w) == w`, so it
# is consistent with the package's `_su2_*` convention with no sign guesswork.
# Rows/cols run m = +s … -s. `twos = 2s`.
function _wignerD(twos::Integer, w::AbstractMatrix)
    n = twos + 1
    T = complex(float(eltype(w)))
    D = zeros(T, n, n)
    ms = twos:-2:-twos
    for (i, mp2) in enumerate(ms), (j, m2) in enumerate(ms)
        A = (twos + mp2) ÷ 2
        B = (twos - mp2) ÷ 2
        C = (twos + m2) ÷ 2
        Dd = (twos - m2) ÷ 2
        off = (mp2 + m2) ÷ 2
        pref = sqrt(factorial(A) * factorial(B) * factorial(C) * factorial(Dd))
        s = zero(T)
        for k in max(0, off):min(A, C)
            denom = factorial(k) * factorial(A - k) * factorial(C - k) * factorial(k - off)
            s += (w[1, 1]^k * w[1, 2]^(A - k) * w[2, 1]^(C - k) * w[2, 2]^(k - off)) / denom
        end
        D[i, j] = pref * s
    end
    return D
end

# Evolve a spin state under one frame step with four-vector map `transform` and its
# explicit 2x2 SL(2,C) representation `U_step`. The Wigner rotation
# w = u(p′)⁻¹ · U_step · u(p) is unitary (rest→rest); boost non-unitarity cancels.
function _evolve_spin(s::SpinState, transform, U_step)
    p′ = transform(s.p)
    u = _prep_su2(s.basis, s.p)
    u′ = _prep_su2(s.basis, p′)
    w = u′ \ (U_step * u)
    c′ = _wignerD(s.twos, w) * s.coeffs
    return SpinState(s.basis, p′, s.twos, c′)
end

"""
    to_basis(s::SpinState, newbasis) -> SpinState

Re-express the same physical state at the same carrier momentum in `newbasis`.
The coefficients are rotated by `Dˢ(u_new(p)⁻¹ · u_old(p))`, the rest-frame rotation
relating the two conventions.
"""
function to_basis(s::SpinState, newbasis::AbstractSpinBasis)
    u_old = _prep_su2(s.basis, s.p)
    u_new = _prep_su2(newbasis, s.p)
    c′ = _wignerD(s.twos, u_new \ u_old) * s.coeffs
    return SpinState(newbasis, s.p, s.twos, c′)
end
