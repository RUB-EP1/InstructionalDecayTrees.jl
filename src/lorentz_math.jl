function _rz_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    z = zero(c)
    o = one(c)
    return [
        c -s z z
        s  c z z
        z z o z
        z z z o
    ]
end

function _ry_xyze(θ::Real)
    c, s = cos(θ), sin(θ)
    z = zero(c)
    o = one(c)
    return [
         c z s z
        z o z z
        -s z c z
        z z z o
    ]
end

function _bz_xyze(ξ::Real)
    g = cosh(ξ)
    bg = sinh(ξ)
    z = zero(g)
    o = one(g)
    return [
        o z z z
        z o z z
        z z g bg
        z z bg g
    ]
end

function _decode_rotation_zyz_xyze(R::AbstractMatrix)
    ϕ = atan(R[2, 3], R[1, 3])
    oneR = one(R[3, 3])
    θ = acos(clamp(R[3, 3], -oneR, oneR))
    ψ = atan(R[3, 2], -R[3, 1])
    return (ϕ = ϕ, θ = θ, ψ = ψ)
end

function _su2_rz(θ::Real)
    h = θ / 2
    Tc = Complex{typeof(float(h))}
    return Tc[
        cis(-h) 0
        0 cis(h)
    ]
end

function _su2_ry(θ::Real)
    h = θ / 2
    c = cos(h)
    s = sin(h)
    Tc = Complex{typeof(float(c))}
    return Tc[
        c -s
        s c
    ]
end

function _su2_bz(ξ::Real)
    h = ξ / 2
    Tc = Complex{typeof(float(h))}
    return Tc[
        exp(h) 0
        0 exp(-h)
    ]
end

_wrap2pi(α::Real) = mod(α + π, 2π) - π

_build_su2(ϕ, θ, ξ, ϕ_rf, θ_rf, ψ_rf) =
    _su2_rz(ϕ) * _su2_ry(θ) * _su2_bz(ξ) *
    _su2_rz(ϕ_rf) * _su2_ry(θ_rf) * _su2_rz(ψ_rf)

normalize_psi(ψ::Real) = mod(ψ + π, 4π) - π

_real_eltype(::Type{T}) where {T<:Real} = T
_real_eltype(::Type{Complex{T}}) where {T<:Real} = T

"""
    _default_atol(x)

Type-scaled tolerance: `10 * eps(T)` for the scalar element type of `x`
(a matrix or its element type).
"""
_default_atol(x) = 10 * eps(_real_eltype(eltype(x)))

function _decode_rotation_zyz_su2(U::AbstractMatrix; atol::Real=_default_atol(U))
    a = U[1, 1]
    c = U[2, 1]
    t = abs(a)
    s = abs(c)
    θ = 2 * atan(s, t)

    if s < atol
        # θ ≈ 0: only (ϕ + ψ) is fixed
        ϕ = zero(real(θ))
        ψ = 2 * angle(U[2, 2])
    elseif t < atol
        # θ ≈ π: only (ϕ - ψ) is fixed
        ϕ = 2 * angle(c)
        ψ = zero(real(θ))
    else
        σ = 2 * angle(U[2, 2]) # ϕ + ψ
        δ = 2 * angle(c)       # ϕ - ψ
        ϕ = (σ + δ) / 2
        ψ = (σ - δ) / 2
    end

    # Rz has period 4π, so wrapping ϕ into (-π, π] can shift it by ±2π, which
    # flips the sign of the rebuilt SU(2) matrix. Shift ψ by the same amount
    # (normalize_psi is 4π-periodic, hence sign-preserving) so that
    # Rz(ϕ)·Ry(θ)·Rz(ψ) reconstructs U exactly, never -U.
    ϕw = _wrap2pi(real(ϕ))
    ψc = real(ψ) + (real(ϕ) - ϕw)
    return (ϕ = ϕw, θ = real(θ), ψ = normalize_psi(ψc))
end

function _decode_boost_xyze(M::AbstractMatrix; atol::Real=_default_atol(M))
    # In (px,py,pz,E) basis, boosting the rest vector [0,0,0,1]
    # corresponds to taking the fourth column.
    v = M[:, 4]
    γ = v[4]
    abs_mom = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    oneγ = one(γ)
    zerom = zero(abs_mom)
    # Products of helicity-frame transforms can be pure rotations up to
    # sqrt(eps) noise in the boost column. Decode those as zero boosts.
    rest_atol = max(atol, sqrt(eps(_real_eltype(eltype(M)))))

    if abs(γ - oneγ) < rest_atol && abs_mom < rest_atol
        zeroγ = zero(γ)
        return (ϕ = zeroγ, θ = zeroγ, ξ = zeroγ)
    end

    if γ < oneγ
        if isapprox(γ, oneγ; atol = atol)
            γ = oneγ
        else
            error("gamma < 1 in Lorentz decode: not a valid Lorentz transformation.")
        end
    end

    ξ = acosh(γ)
    ϕ = atan(v[2], v[1])
    cinput = abs(abs_mom) <= atol ? zerom : v[3] / abs_mom
    θ = acos(clamp(cinput, -oneγ, oneγ))

    if abs(γ - oneγ) < atol
        zeroγ = zero(γ)
        return (ϕ = zeroγ, θ = zeroγ, ξ = zeroγ)
    end
    return (ϕ = ϕ, θ = θ, ξ = ξ)
end

function _decode_lorentz_helicity_zyz_xyze(M::AbstractMatrix; atol::Real=_default_atol(M))
    b = _decode_boost_xyze(M; atol = atol)
    M_rf = _bz_xyze(-b.ξ) * _ry_xyze(-b.θ) * _rz_xyze(-b.ϕ) * M
    rot = abs(b.ξ) < atol ? _decode_rotation_zyz_xyze(M[1:3, 1:3]) :
          _decode_rotation_zyz_xyze(M_rf[1:3, 1:3])
    return (ϕ = b.ϕ, θ = b.θ, ξ = b.ξ, ϕ_rf = rot.ϕ, θ_rf = rot.θ, ψ_rf = rot.ψ)
end
