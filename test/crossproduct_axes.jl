using InstructionalDecayTrees
using FourVectors
using LinearAlgebra
using Random
using Test

const IDT = InstructionalDecayTrees

comps(p::FourVector) = [p.px, p.py, p.pz, p.E]
wrap_delta(a, b) = mod(a - b + π, 2π) - π

# A generic four-body configuration in the rest frame of the total, generated
# through a D* / X cascade so the intermediate systems are physical.
function random_objs(rng)
    two_body(M, m1, m2) = begin
        p = sqrt((M^2 - (m1 + m2)^2) * (M^2 - (m1 - m2)^2)) / (2M)
        cosθ = 2rand(rng) - 1
        sinθ = sqrt(1 - cosθ^2)
        ϕ = 2π * rand(rng)
        v = (p * sinθ * cos(ϕ), p * sinθ * sin(ϕ), p * cosθ)
        (FourVector(v...; E = sqrt(p^2 + m1^2)),
            FourVector((-).(v)...; E = sqrt(p^2 + m2^2)))
    end
    from_rest(q, P) = boost_to_rest(q, FourVector(-P.px, -P.py, -P.pz; E = P.E))

    M, m1, m2, m3, m4 = 5.28, 1.86, 0.14, 1.87, 0.49
    mDst = m1 + m2 + rand(rng) * (M - m4 - m3 - m1 - m2)
    mX = mDst + m3 + rand(rng) * (M - m4 - mDst - m3)
    X, p4 = two_body(M, mX, m4)
    Dst_x, D_x = two_body(mX, mDst, m3)
    Dst, p3 = from_rest(Dst_x, X), from_rest(D_x, X)
    a_x, b_x = two_body(mDst, m1, m2)
    p1, p2 = from_rest(a_x, Dst), from_rest(b_x, Dst)
    # random overall orientation, one rotation for the whole event
    α, β, γ = 2π * rand(rng), acos(2rand(rng) - 1), 2π * rand(rng)
    R = p -> p |> Rz(α) |> Ry(β) |> Rz(γ)
    return map(R, (p1, p2, p3, p4))
end

@testset "boost_to_rest is a pure boost" begin
    rng = MersenneTwister(1)
    objs = random_objs(rng)
    P = objs[1] + objs[2] + objs[3]
    θ, ϕ, γ = polar_angle(P), azimuthal_angle(P), boost_gamma(P)

    # the sandwich, spelled out
    for p in objs
        expected = p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ) |> Ry(θ) |> Rz(ϕ)
        @test comps(boost_to_rest(p, P)) ≈ comps(expected) atol = 1e-14
    end

    # it brings P to rest and conserves every mass
    P_rest = boost_to_rest(P, P)
    @test hypot(P_rest.px, P_rest.py, P_rest.pz) < 1e-12
    for p in objs
        @test mass(boost_to_rest(p, P)) ≈ mass(p) atol = 1e-12
    end
end

@testset "ToHelicityFrame == ToRestFrame then realignment" begin
    rng = MersenneTwister(2)
    objs = random_objs(rng)
    P = objs[1] + objs[2]
    θ, ϕ = polar_angle(P), azimuthal_angle(P)
    for p in objs
        realigned = boost_to_rest(p, P) |> Rz(-ϕ) |> Ry(-θ)
        @test comps(transform_to_cmf(p, P)) ≈ comps(realigned) atol = 1e-14
    end
end

@testset "ToRestFrame is the identity for a system already at rest" begin
    objs = random_objs(MersenneTwister(3))
    P = sum(objs)                       # total, at rest by construction
    (moved, _) = apply_decay_instruction((ToRestFrame((1, 2, 3, 4)),), objs)
    for (a, b) in zip(objs, moved)
        @test comps(a) ≈ comps(b) atol = 1e-14
    end
end

@testset "axis markers are inert under a pure boost" begin
    objs = random_objs(MersenneTwister(4))
    state = with_helicity_axes(objs)
    (state, _) = apply_decay_instruction((PlantLabAxes(5),), state)
    planted = helicity_axes_at(state, 5)
    (state, _) = apply_decay_instruction((ToRestFrame((1, 2)),), state)
    @test helicity_axes_at(state, 5).ẑ === planted.ẑ
    @test helicity_axes_at(state, 5).x̂ === planted.x̂
end

@testset "euler_zxz reduces to (ϕ, θ) against the coordinate axes" begin
    rng = MersenneTwister(5)
    for _ in 1:200
        v = FourVector(randn(rng), randn(rng), randn(rng); E = 10.0)
        ang = euler_zxz((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (v.px, v.py, v.pz))
        @test abs(wrap_delta(ang.α, azimuthal_angle(v))) < 1e-12
        @test ang.β ≈ polar_angle(v) atol = 1e-12
        @test ang.γ == 0.0
    end
end

@testset "carried axes vs realigned frame: same angles" begin
    rng = MersenneTwister(6)
    # (((1,2),3),4) and ((1,2),(3,4)) topologies, every vertex.
    walks = [
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)),
            MeasureEulerZXZ(:a, (1, 2, 3), 5)),
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), TransportAxes(5, (1, 2, 3)),
            ToRestFrame((1, 2, 3)), MeasureEulerZXZ(:a, (1, 2), 5)),
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), TransportAxes(5, (1, 2, 3)),
            ToRestFrame((1, 2, 3)), TransportAxes(5, (1, 2)), ToRestFrame((1, 2)),
            MeasureEulerZXZ(:a, (1,), 5)),
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), MeasureEulerZXZ(:a, (1, 2), 5)),
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), TransportAxes(5, (1, 2)),
            ToRestFrame((1, 2)), MeasureEulerZXZ(:a, (1,), 5)),
        # second daughter of the root
        (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), TransportAxes(5, (3, 4)),
            ToRestFrame((3, 4)), MeasureEulerZXZ(:a, (3,), 5)),
    ]
    # The events are already in the total rest frame, so the realigning programs
    # start in the current frame (no leading ToHelicityFrame at the root).
    realigned = [
        (MeasureCosThetaPhi(:a, (1, 2, 3)),),
        (ToHelicityFrame((1, 2, 3)), MeasureCosThetaPhi(:a, (1, 2))),
        (ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2)), MeasureCosThetaPhi(:a, (1,))),
        (MeasureCosThetaPhi(:a, (1, 2)),),
        (ToHelicityFrame((1, 2)), MeasureCosThetaPhi(:a, (1,))),
        (ToHelicityFrame((3, 4)), MeasureCosThetaPhi(:a, (3,))),
    ]

    worst = 0.0
    for _ in 1:200
        objs = random_objs(rng)
        for (walk, plain) in zip(walks, realigned)
            (_, rw) = apply_decay_instruction(walk, with_helicity_axes(objs))
            (_, rp) = apply_decay_instruction(plain, objs)
            worst = max(worst,
                abs(wrap_delta(rw.a.α, rp.a.ϕ)),
                abs(rw.a.β - acos(clamp(rp.a.cosθ, -1, 1))))
        end
    end
    @test worst < 1e-10
end

@testset "measurement result schema and branch folding" begin
    objs = random_objs(MersenneTwister(7))
    program = (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)), MeasureEulerZXZ(:v, (1, 2), 5))
    (_, res) = apply_decay_instruction(program, with_helicity_axes(objs))
    @test keys(res.v) == (:α, :β, :γ, :cosβ)
    @test res.v.cosβ ≈ cos(res.v.β)
    @test -π <= res.v.α < π

    program2 = (PlantLabAxes(5), ToRestFrame((1, 2, 3, 4)),
        MeasureEulerZXZ(:v, (1, 2), 5; branch = 2))
    (_, res2) = apply_decay_instruction(program2, with_helicity_axes(objs))
    @test -2π <= res2.v.α < 0
    @test abs(wrap_delta(res2.v.α, res.v.α)) < 1e-12   # same angle, other branch
    @test res2.v.β ≈ res.v.β
end

@testset "axis slots are rejected as momenta" begin
    objs = with_helicity_axes(random_objs(MersenneTwister(8)))
    @test_throws ArgumentError apply_decay_instruction((ToRestFrame((1, 5)),), objs)
end

@testset "tracking through ToRestFrame" begin
    objs = random_objs(MersenneTwister(9))
    state = init_tracked_state(with_helicity_axes(objs))
    program = (PlantLabAxes(5), ToRestFrame((1, 2, 3)), TransportAxes(5, (1, 2)))
    (final, _) = apply_decay_instruction(program, state)

    # the tracked matrix must reproduce the transformation actually applied
    for k in 1:4
        col = final.tracker.Λ * [objs[k].px, objs[k].py, objs[k].pz, objs[k].E]
        @test col ≈ [final.objs[k].px, final.objs[k].py, final.objs[k].pz, final.objs[k].E] atol = 1e-12
    end

    # a pure boost and the corresponding helicity frame differ by a rotation only
    cmp = compare_instruction_paths(
        (ToRestFrame((1, 2, 3)),), (ToHelicityFrame((1, 2, 3)),), objs)
    @test abs(decode_lorentz_helicity(cmp.relative).ξ) < 1e-10
end

@testset "compact display" begin
    @test repr(PlantLabAxes(5)) == "PlantLabAxes(5)"
    @test repr(ToRestFrame((1, 2))) == "ToRestFrame((1, 2))"
    @test repr(ToRestFrame(1)) == "ToRestFrame(1)"
    @test repr(TransportAxes(5, (1, 2))) == "TransportAxes(5, (1, 2))"
    @test repr(TransportAxes(5, 1)) == "TransportAxes(5, 1)"
    @test repr(MeasureEulerZXZ(:v1, (1, 2), 5)) == "MeasureEulerZXZ(:v1, (1, 2), 5)"
    @test repr(MeasureEulerZXZ(:v1, 1, 5; branch = 2)) ==
          "MeasureEulerZXZ(:v1, 1, 5; branch=2)"
end
