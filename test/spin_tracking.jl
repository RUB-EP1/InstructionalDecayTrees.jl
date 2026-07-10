using InstructionalDecayTrees
using FourVectors
using LinearAlgebra
using Test

const IDT = InstructionalDecayTrees

# Shared fixtures
const P1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
const P2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
const P3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
const P4 = FourVector(0.25, -0.1, 0.35; M = 0.35)

# Wigner rotation acting on a spin state's coefficients over one step, computed
# directly from the accumulated SU(2) `U` and the basis preparations.
_wig(basis, twos, p0, pf, U) = IDT._wignerD(twos, IDT._prep_su2(basis, pf) \ (U * IDT._prep_su2(basis, p0)))

@testset "_wignerD representation" begin
    w = IDT._su2_rz(0.3) * IDT._su2_ry(0.7) * IDT._su2_bz(0.2)   # generic SL(2,C)
    w2 = IDT._su2_ry(1.1) * IDT._su2_rz(-0.4)

    @test IDT._wignerD(1, w) ≈ w atol = 1e-14   # spin-1/2 rep is the matrix itself

    for twos in (1, 2, 3, 4)
        @test IDT._wignerD(twos, Matrix{ComplexF64}(I, 2, 2)) ≈
              Matrix{ComplexF64}(I, twos + 1, twos + 1) atol = 1e-14
        @test IDT._wignerD(twos, w * w2) ≈ IDT._wignerD(twos, w) * IDT._wignerD(twos, w2) atol =
            1e-12
        u = IDT._su2_rz(0.9) * IDT._su2_ry(2.1) * IDT._su2_rz(-1.3)   # unitary -> unitary
        D = IDT._wignerD(twos, u)
        @test D * D' ≈ Matrix{ComplexF64}(I, twos + 1, twos + 1) atol = 1e-12
    end

    α = 0.5
    Drz = IDT._wignerD(2, IDT._su2_rz(α))
    @test Drz ≈ Diagonal([cis(-1 * α), cis(0 * α), cis(1 * α)]) atol = 1e-14

    β = 0.7   # spin-1 small-d against analytic Wigner d^1(β)
    c, s = cos(β), sin(β)
    d1 = [
        (1+c)/2 -s/sqrt(2) (1-c)/2
        s/sqrt(2) c -s/sqrt(2)
        (1-c)/2 s/sqrt(2) (1+c)/2
    ]
    @test IDT._wignerD(2, IDT._su2_ry(β)) ≈ d1 atol = 1e-12
end

@testset "PlaneAlign explicit SU(2) matches its SO(3)" begin
    axis_z = FourVector(0.3, -0.4, 0.5; E = 2.0)
    axis_x = FourVector(-0.6, 0.2, 0.1; E = 1.5)
    transform = p -> rotate_to_plane(p, axis_z, axis_x)
    M = IDT._step_matrix(transform, Float64)

    ϕ_z = azimuthal_angle(axis_z)
    θ_z = polar_angle(axis_z)
    α = azimuthal_angle(axis_x |> Rz(-ϕ_z) |> Ry(-θ_z))
    R = IDT._rz_xyze(-α) * IDT._ry_xyze(-θ_z) * IDT._rz_xyze(-ϕ_z)
    @test M ≈ R atol = 1e-12
end

@testset "Spin state basics" begin
    s = spin_state(Helicity(), P1, 1, 1 // 2)
    @test s.twos == 1 && s.coeffs == ComplexF64[1, 0]
    @test spin_state(Canonical(), P1, 2, -1).coeffs == ComplexF64[0, 0, 1]
    @test projections(1) == (1 // 2, -1 // 2)
    @test projections(2) == (1 // 1, 0 // 1, -1 // 1)
    @test_throws ArgumentError spin_state(Helicity(), P1, 1, 1)  # 2m parity mismatch
    @test_throws ArgumentError spin_state(Helicity(), P1, 1, 3 // 2)  # |m| > s
end

# ---------------------------------------------------------------------------
# Helicity basis: physical properties
# ---------------------------------------------------------------------------
@testset "Helicity: collinear boost is invariant" begin
    # Carrier and boosted system both along +z: no Wigner rotation, coeffs fixed.
    q1 = FourVector(0.0, 0.0, 0.3; M = 0.2)
    q2 = FourVector(0.0, 0.0, 0.5; M = 0.4)
    c0 = normalize(ComplexF64[0.6, 0.8])
    st = init_tracked_state((q1, q2); spins = (a = spin_state(Helicity(), q1, c0),))
    (fin, _) = apply_decay_instruction((ToHelicityFrame((1, 2)),), st)
    @test fin.spins.a.coeffs ≈ c0 atol = 1e-12
end

@testset "Helicity: pure rotation conserves helicity (diagonal Wigner rot.)" begin
    # Under a pure rotation the massive-helicity little group is Rz(γ), so the
    # Wigner rotation is diagonal: helicity magnitudes are preserved per component.
    objs = (P1, P2, P3)
    for twos in (1, 2)
        c0 = normalize(ComplexF64[1.0 + 0.2im, -0.3, 0.5im][1:(twos+1)])
        st = init_tracked_state(objs; spins = (a = spin_state(Helicity(), P3, c0),))
        (fin, _) = apply_decay_instruction((PlaneAlign(1, 2),), st)

        w = IDT._prep_su2(Helicity(), fin.objs[3]) \ (fin.tracker.U * IDT._prep_su2(Helicity(), P3))
        @test w ≈ Diagonal(diag(w)) atol = 1e-12               # w is diagonal
        @test abs.(fin.spins.a.coeffs) ≈ abs.(c0) atol = 1e-11 # magnitudes preserved
        @test norm(fin.spins.a.coeffs) ≈ 1 atol = 1e-11
    end
end

# ---------------------------------------------------------------------------
# Canonical basis: physical properties
# ---------------------------------------------------------------------------
@testset "Canonical: pure rotation acts as D(R), momentum-independent" begin
    # Defining property of the canonical (rotationless) boost: under a pure
    # rotation R the spin follows by exactly D(R), regardless of the carrier
    # momentum. Verify identical D(R)·c0 for carriers on two different particles.
    objs = (P1, P2, P3)
    for twos in (1, 2)
        c0 = normalize(ComplexF64[0.5 + 0.1im, 0.8, -0.2][1:(twos+1)])
        results = map((1, 3)) do carrier
            st = init_tracked_state(objs; spins = (a = spin_state(Canonical(), objs[carrier], c0),))
            (fin, _) = apply_decay_instruction((PlaneAlign(1, 2),), st)
            (fin.spins.a.coeffs, fin.tracker.U)
        end
        DR = IDT._wignerD(twos, results[1][2])   # D of the pure rotation U_step
        @test results[1][1] ≈ DR * c0 atol = 1e-11   # carrier 1
        @test results[2][1] ≈ DR * c0 atol = 1e-11   # carrier 3 (same D(R))
    end
end

@testset "Canonical: non-collinear boost induces mixing (Wigner rotation)" begin
    # A non-collinear boost rotates a canonical eigenstate into a superposition.
    q1 = FourVector(0.5, 0.0, 0.0; M = 0.2)   # along +x
    q2 = FourVector(0.0, 0.0, 0.6; M = 0.4)   # boost axis along +z (non-collinear)
    st = init_tracked_state((q1, q2); spins = (a = spin_state(Canonical(), q1, 1, 1 // 2),))
    (fin, _) = apply_decay_instruction((ToHelicityFrame((1, 2)),), st)
    @test abs(fin.spins.a.coeffs[2]) > 1e-3          # lower component populated
    @test norm(fin.spins.a.coeffs) ≈ 1 atol = 1e-12  # but still unitary
end

# ---------------------------------------------------------------------------
# Helicity <-> Canonical relations
# ---------------------------------------------------------------------------
@testset "Relation: to_basis(H->C) = D(R(ϕ,θ))" begin
    # At fixed p, u_C⁻¹ u_H = R(ϕ,θ), so canonical coeffs = D(R(ϕ,θ)) · helicity.
    for twos in (1, 2, 3)
        c0 = normalize(ComplexF64[0.7, -0.4 + 0.3im, 0.5, -0.2im][1:(twos+1)])
        sh = spin_state(Helicity(), P2, c0)
        sc = to_basis(sh, Canonical())
        ϕ, θ = azimuthal_angle(P2), polar_angle(P2)
        R = IDT._su2_rz(ϕ) * IDT._su2_ry(θ)
        @test sc.coeffs ≈ IDT._wignerD(twos, R) * c0 atol = 1e-12
        @test to_basis(sc, Helicity()).coeffs ≈ c0 atol = 1e-12   # roundtrip
    end
end

@testset "Relation: bases coincide for momentum along +z" begin
    pz = FourVector(0.0, 0.0, 0.5; M = 0.3)   # θ = 0, R(ϕ,θ) = I
    c0 = normalize(ComplexF64[0.6, 0.8])
    @test to_basis(spin_state(Helicity(), pz, c0), Canonical()).coeffs ≈ c0 atol = 1e-12
end

@testset "Relation: evolve-in-each-basis then convert agree" begin
    # Same physical state (linked by to_basis), evolved through the same path,
    # remains the same physical state.
    objs = (P1, P2, P3)
    path = (ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2)))
    c0 = normalize(ComplexF64[0.7, -0.4 + 0.3im])
    sh = spin_state(Helicity(), P1, c0)
    sc = to_basis(sh, Canonical())

    (fin_h, _) = apply_decay_instruction(path, init_tracked_state(objs; spins = (a = sh,)))
    (fin_c, _) = apply_decay_instruction(path, init_tracked_state(objs; spins = (a = sc,)))
    @test to_basis(fin_c.spins.a, Helicity()).coeffs ≈ fin_h.spins.a.coeffs atol = 1e-11
end

# ---------------------------------------------------------------------------
# Consistency with the accumulated Lorentz tracker (both bases)
# ---------------------------------------------------------------------------
@testset "Tracker consistency (telescoping to accumulated U)" begin
    objs = (P1, P2, P3)
    path = (ToHelicityFrame((1, 2, 3)), PlaneAlign(2, 3), ToHelicityFrame((2, 3)))
    for basis in (Helicity(), Canonical()), twos in (1, 2)
        c0 = normalize(ComplexF64[1.0 + 0.2im, -0.3, 0.5im][1:(twos+1)])
        st = init_tracked_state(objs; spins = (a = spin_state(basis, P1, c0),))
        (fin, _) = apply_decay_instruction(path, st)
        expected = _wig(basis, twos, P1, fin.objs[1], fin.tracker.U) * c0
        @test fin.spins.a.coeffs ≈ expected atol = 1e-10
        @test norm(fin.spins.a.coeffs) ≈ 1 atol = 1e-10
    end
end

@testset "Spin tracking through Gottfried-Jackson frame" begin
    objs = (P1, P2, P3, P4)
    c0 = normalize(ComplexF64[0.5 + 0.1im, 0.8])
    st = init_tracked_state(objs; spins = (a = spin_state(Helicity(), P1, c0),))
    (fin, _) = apply_decay_instruction((ToGottfriedJacksonFrame((2, 3, 4), 2, 4),), st)
    @test norm(fin.spins.a.coeffs) ≈ 1 atol = 1e-10
    @test fin.spins.a.coeffs ≈ _wig(Helicity(), 1, P1, fin.objs[1], fin.tracker.U) * c0 atol = 1e-10
end
