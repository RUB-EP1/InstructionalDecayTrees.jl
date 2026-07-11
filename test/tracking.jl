using LinearAlgebra

@testset "TrackedState dispatch" begin
    p1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
    p2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
    p3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
    objs_local = (p1, p2, p3)

    program = (
        ToHelicityFrame((1, 2, 3)),
        MeasureInvariant(:m2, (1, 2)),
        ToHelicityFrame((1, 2)),
        MeasurePolar(:theta1, 1),
    )

    (final_plain, res_plain) = apply_decay_instruction(program, objs_local)
    (final_tracked, res_tracked) = apply_decay_instruction(program, init_tracked_state(objs_local))

    @test final_tracked isa TrackedState
    @test res_plain == res_tracked

    for i in eachindex(final_plain)
        @test final_plain[i].E ≈ final_tracked.objs[i].E atol = 1e-12
        @test final_plain[i].px ≈ final_tracked.objs[i].px atol = 1e-12
        @test final_plain[i].py ≈ final_tracked.objs[i].py atol = 1e-12
        @test final_plain[i].pz ≈ final_tracked.objs[i].pz atol = 1e-12
    end
end

@testset "Instruction path comparison" begin
    p1 = FourVector(0.5, 0.1, -0.2; M = 0.14)
    p2 = FourVector(-0.4, 0.2, 0.3; M = 0.30)
    p3 = FourVector(-0.1, -0.3, -0.1; M = 0.20)
    objs_local = (p1, p2, p3)

    path = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrame((1, 2)),
    )

    cmp = compare_instruction_paths(path, path, objs_local)
    I4 = Matrix{Float64}(I, 4, 4)

    @test cmp.tracker1.Λ ≈ cmp.tracker2.Λ atol = 1e-12
    @test cmp.relative.Λ ≈ I4 atol = 1e-10
end

@testset "TrackedState Float32 support" begin
    p1 = FourVector(Float32(0.3), Float32(0.2), Float32(0.1); M = Float32(0.2))
    p2 = FourVector(Float32(-0.1), Float32(0.4), Float32(-0.3); M = Float32(0.4))
    p3 = FourVector(Float32(-0.2), Float32(-0.6), Float32(0.2); M = Float32(0.3))
    objs_local = (p1, p2, p3)

    path = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrame((1, 2)),
    )

    (state, _) = apply_decay_instruction(path, init_tracked_state(objs_local; T = Float32))
    @test state isa TrackedState
    @test eltype(state.tracker.Λ) == Float32

    cmp = compare_instruction_paths(path, path, objs_local; T = Float32)
    atol32 = InstructionalDecayTrees._default_atol(cmp.relative.Λ)
    @test atol32 ≈ 10 * eps(Float32)
    @test cmp.relative.Λ ≈ Matrix{Float32}(I, 4, 4) atol = atol32
    @test wigner_zyz(cmp.relative) isa NamedTuple
    @test InstructionalDecayTrees._wigner_zyz_su2(cmp.relative) isa NamedTuple
end

@testset "Explicit SU2 step uses rapidity (not gamma)" begin
    p1 = FourVector(0.25, -0.15, 0.30; M = 0.2)
    p2 = FourVector(-0.10, 0.05, -0.20; M = 0.3)
    p3 = FourVector(0.15, 0.10, -0.05; M = 0.25)
    objs_local = (p1, p2, p3)
    st = init_tracked_state(objs_local)

    P_tot = p1 + p2 + p3
    ϕ = azimuthal_angle(P_tot)
    θ = polar_angle(P_tot)
    γ = boost_gamma(P_tot)
    ξ = acosh(γ)

    transform = p -> transform_to_cmf(p, P_tot)
    U_inferred = begin
        M_step = InstructionalDecayTrees._step_matrix(transform, Float64)
        d = InstructionalDecayTrees._decode_lorentz_helicity_zyz_xyze(M_step)
        InstructionalDecayTrees._build_su2(d.ϕ, d.θ, d.ξ, d.ϕ_rf, d.θ_rf, d.ψ_rf)
    end
    U_explicit = InstructionalDecayTrees._su2_bz(-ξ) * InstructionalDecayTrees._su2_ry(-θ) *
                 InstructionalDecayTrees._su2_rz(-ϕ)

    err_plus = sum(abs2, U_explicit .- U_inferred)
    err_minus = sum(abs2, U_explicit .+ U_inferred)
    @test min(err_plus, err_minus) < 1e-10
end

@testset "Cross-check against decayangle (4-body, particle 3)" begin
    # Same four-vectors used in existing tests and Python cross-check harness.
    p1 = FourVector(-0.1467, 0.2235, -0.7847; E = 2.0452)
    p2 = FourVector(-0.0873, 0.1803, -0.5584; E = 0.7718)
    p3 = FourVector(0.0056, -0.0349, 0.1413; E = 0.2017)
    p4 = FourVector(0.2284, -0.3689, 1.2019; E = 2.2606)
    objs_local = (p1, p2, p3, p4)

    # Reference topology ((12)3)4 path to particle 3:
    path_ref = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrameParticle2(3),
    )

    # Other topology 1(2(34)) path to particle 3:
    path_other = (
        ToHelicityFrameParticle2((2, 3, 4)),
        ToHelicityFrameParticle2((3, 4)),
        ToHelicityFrame((3,)),
    )

    cmp = compare_instruction_paths(path_ref, path_other, objs_local)

    # Expected relative matrix from decayangle:
    # rel = boost_other @ inv(boost_reference), in (px,py,pz,E) basis.
    rel_expected = [
        0.9262423551164856 -0.023396446110968982 -0.3762016824758271 0.0
        0.02166279138320244 0.9997262650053029 -0.008838468647614112 0.0
        0.37630549166568683 3.6985450305581945e-5 0.9264956425015303 0.0
        0.0 0.0 0.0 1.0
    ]

    @test cmp.relative.Λ ≈ rel_expected atol = 5e-10

    ang = wigner_zyz(cmp.relative)
    @test ang.ϕ ≈ -3.1181030111491173 atol = 5e-10
    @test ang.θ ≈ 0.38580543795548133 atol = 5e-10
    @test ang.ψ ≈ 9.42467967506533 atol = 5e-10
end

@testset "SU2 branch resolves 2π phase" begin
    # Non-singular ZYZ rotation with explicit spinor sign flip.
    # Same Λ, opposite U branch should force ψ to the shifted 2π branch.
    ϕ0 = 0.4
    θ0 = 1.1
    ψ0 = 0.7
    Λ = InstructionalDecayTrees._rz_xyze(ϕ0) *
        InstructionalDecayTrees._ry_xyze(θ0) *
        InstructionalDecayTrees._rz_xyze(ψ0)
    U = -InstructionalDecayTrees._build_su2(0.0, 0.0, 0.0, ϕ0, θ0, ψ0)
    t = LorentzTracker(Λ, U)

    ang = wigner_zyz(t)
    @test ang.ϕ ≈ ϕ0 atol = 1e-12
    @test ang.θ ≈ θ0 atol = 1e-12
    @test ang.ψ ≈ (ψ0 + 2π) atol = 1e-12
end

@testset "PlaneAlign tracks an explicit phase-correct SU2" begin
    # Regression for issue #27: the inferred-U_step fallback rebuilds the SU(2)
    # from the SO(3) block, whose ±1 branch is an atan/acos artifact and flips
    # discontinuously under smooth kinematics.
    su2_explicit(axis_z, axis_x) = begin
        ϕ_z = azimuthal_angle(axis_z)
        θ_z = polar_angle(axis_z)
        α = azimuthal_angle(axis_x |> Rz(-ϕ_z) |> Ry(-θ_z))
        InstructionalDecayTrees._su2_rz(-α) *
        InstructionalDecayTrees._su2_ry(-θ_z) *
        InstructionalDecayTrees._su2_rz(-ϕ_z)
    end
    tracked_U(objs) = begin
        (st, _) = apply_decay_instruction((PlaneAlign(1, 2),), init_tracked_state(objs))
        st.tracker
    end

    q2 = FourVector(-0.6, 0.2, 0.1; E = 1.5)
    q3 = FourVector(0.1, 0.3, -0.2; E = 1.0)

    # continuity: before the fix the tracked U flipped sign across this ε-step
    # (ϕ_z ≈ 2.8198…) while Λ stayed continuous
    for ϕ in (2.819841, 2.819843)
        q1 = FourVector(0.5 * cos(ϕ), 0.5 * sin(ϕ), 0.4; E = 2.0)
        t = tracked_U((q1, q2, q3))
        @test t.U ≈ su2_explicit(q1, q2) atol = 1e-10   # exact, no ± ambiguity
    end

    # generic configuration: exactness and SO(3)-image consistency with Λ
    q1 = FourVector(0.3, -0.4, 0.5; E = 2.0)
    t = tracked_U((q1, q2, q3))
    @test t.U ≈ su2_explicit(q1, q2) atol = 1e-12
    σs = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
    R_adj = [real(0.5 * tr(σi * t.U * σj * t.U')) for σi in σs, σj in σs]
    @test R_adj ≈ t.Λ[1:3, 1:3] atol = 1e-10
end

@testset "SO3 and SU2 Wigner decoders agree on pure rotations" begin
    wrap2(a, b) = mod(a - b + π, 2π) - π
    wrap4(a, b) = mod(a - b + 2π, 4π) - 2π
    for (ϕ0, θ0, ψ0) in ((0.4, 1.1, 0.7), (-2.3, 0.9, 2.6), (1.9, 2.2, -1.4))
        Λ = InstructionalDecayTrees._rz_xyze(ϕ0) *
            InstructionalDecayTrees._ry_xyze(θ0) *
            InstructionalDecayTrees._rz_xyze(ψ0)
        U = InstructionalDecayTrees._build_su2(0.0, 0.0, 0.0, ϕ0, θ0, ψ0)
        t = LorentzTracker(Λ, U)

        a_so3 = InstructionalDecayTrees._wigner_zyz_so3(t)
        a_su2 = InstructionalDecayTrees._wigner_zyz_su2(t)
        a_mix = wigner_zyz(t)

        @test abs(wrap2(a_so3.ϕ, a_su2.ϕ)) < 1e-12
        @test abs(wrap2(a_so3.θ, a_su2.θ)) < 1e-12
        @test abs(wrap4(a_so3.ψ, a_su2.ψ)) < 1e-12

        @test abs(wrap2(a_mix.ϕ, a_su2.ϕ)) < 1e-12
        @test abs(wrap2(a_mix.θ, a_su2.θ)) < 1e-12
        @test abs(wrap4(a_mix.ψ, a_su2.ψ)) < 1e-12
    end
end

@testset "Public Wigner API handles near-pure small-mass rotations" begin
    wrap2(a) = mod(a + π, 2π) - π

    q1(m) = FourVector(0.20, -0.10, -0.50; M = m)
    q2 = FourVector(0.35, 0.22, 0.28; M = 0.45)
    q3 = FourVector(-0.18, 0.31, -0.12; M = 0.60)
    q4 = FourVector(-0.37, -0.43, 0.34; M = 0.80)

    path_123_12_1 = (
        ToHelicityFrame((1, 2, 3)),
        ToHelicityFrame((1, 2)),
        ToHelicityFrame(1),
    )
    path_13_1 = (ToHelicityFrame((1, 3)), ToHelicityFrame(1))

    cmp = compare_instruction_paths(path_123_12_1, path_13_1, (q1(1e-3), q2, q3, q4))
    ang = wigner_zyz(cmp.relative)

    @test ang.θ / π ≈ 0.0004819884356761267 atol = 1e-8
    @test abs(wrap2(ang.ϕ + ang.ψ)) / π ≈ 0.7685257868157153 atol = 1e-8
end
