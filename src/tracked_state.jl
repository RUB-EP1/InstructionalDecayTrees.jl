"""
    TrackedState(objs, tracker::LorentzTracker)

Carries both physics objects and an accumulated Lorentz tracker through
`apply_decay_instruction` dispatch.
"""
struct TrackedState{O,L}
    objs::O
    tracker::L
end

"""
    init_tracked_state(objs; T=Float64)

Create a [`TrackedState`](@ref) with `objs` and an identity [`LorentzTracker`](@ref).

Use this when you want `apply_decay_instruction` to return both transformed
objects and the accumulated Lorentz transformation.
"""
init_tracked_state(objs; T::Type{<:Real}=Float64) = TrackedState(objs, LorentzTracker(T))

_as_column(p::FourVector) = [p.px, p.py, p.pz, p.E]

function _basis4(::Type{T}) where {T<:Real}
    return (
        FourVector(one(T), zero(T), zero(T); E = zero(T)),
        FourVector(zero(T), one(T), zero(T); E = zero(T)),
        FourVector(zero(T), zero(T), one(T); E = zero(T)),
        FourVector(zero(T), zero(T), zero(T); E = one(T)),
    )
end

function _step_matrix(transform, ::Type{T}) where {T<:Real}
    basis = _basis4(T)
    M = Matrix{T}(undef, 4, 4)
    for (j, b) in enumerate(basis)
        bp = transform(b)
        M[:, j] = _as_column(bp)
    end
    return M
end

function _apply_step_with_tracking(state::TrackedState, transform; U_step)
    new_objs = map(transform, state.objs)
    Tobj = typeof(first(state.objs).px)
    # The 4x4 side is generic: the linear map is inferred by transforming basis
    # vectors. The SU(2) side cannot be inferred — SU(2)→SO(3) is 2:1, so the
    # vector action does not determine the spinor sign — hence U_step is a
    # required keyword: every frame instruction must construct its SU(2)
    # explicitly from the physical angles (issue #27).
    M_step = _step_matrix(transform, Tobj)
    new_tracker = LorentzTracker(M_step * state.tracker.Λ, U_step * state.tracker.U)
    return TrackedState(new_objs, new_tracker)
end

function apply_decay_instruction(instr::ToHelicityFrame, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)
    ϕ = azimuthal_angle(P_tot)
    θ = polar_angle(P_tot)
    γ = boost_gamma(P_tot)
    # FourVectors uses gamma-factor parameterization for Bz; SU2 boost uses rapidity.
    ξ = acosh(γ)
    transform = p -> transform_to_cmf(p, P_tot)
    U_step = _su2_bz(-ξ) * _su2_ry(-θ) * _su2_rz(-ϕ)
    return (_apply_step_with_tracking(state, transform; U_step = U_step), (;))
end

function apply_decay_instruction(instr::ToHelicityFrameParticle2, state::TrackedState)
    P_tot = get_fourvector(state.objs, instr.indices)

    P_inv = FourVector(-P_tot.px, -P_tot.py, -P_tot.pz; E = P_tot.E)
    ϕ_inv = azimuthal_angle(P_inv)
    θ_inv = polar_angle(P_inv)
    γ = boost_gamma(P_tot)
    # FourVectors uses gamma-factor parameterization for Bz; SU2 boost uses rapidity.
    ξ = acosh(γ)

    transform = p -> p |> Rz(-ϕ_inv) |> Ry(-θ_inv) |> Ry(-π) |> Bz(-γ)
    U_step = _su2_bz(-ξ) * _su2_ry(-π) * _su2_ry(-θ_inv) * _su2_rz(-ϕ_inv)
    return (_apply_step_with_tracking(state, transform; U_step), (;))
end

function apply_decay_instruction(instr::PlaneAlign, state::TrackedState)
    axis_z = get_fourvector(state.objs, instr.z_idx)
    axis_x = get_fourvector(state.objs, instr.x_idx)
    transform = p -> rotate_to_plane(p, axis_z, axis_x)
    # Explicit, phase-correct SU(2) for the plane-alignment rotation. The angles
    # come from the axis vectors themselves, so the spinor ±1 branch is pinned
    # and varies continuously with the kinematics; the inferred-U_step fallback
    # would rebuild it from the SO(3) block, where the branch is an artifact of
    # atan/acos ranges (issue #27).
    ϕ_z = azimuthal_angle(axis_z)
    θ_z = polar_angle(axis_z)
    α = azimuthal_angle(axis_x |> Rz(-ϕ_z) |> Ry(-θ_z))
    U_step = _su2_rz(-α) * _su2_ry(-θ_z) * _su2_rz(-ϕ_z)
    return (_apply_step_with_tracking(state, transform; U_step = U_step), (;))
end

function apply_decay_instruction(instr::ToGottfriedJacksonFrame, state::TrackedState)
    hel_instr = ToHelicityFrame(instr.system_indices)
    (after_boost, _) = apply_decay_instruction(hel_instr, state)

    target_idx_negated = Tuple(-idx for idx in instr.target_idx)
    plane_instr = PlaneAlign(instr.beam_idx, target_idx_negated)
    return apply_decay_instruction(plane_instr, after_boost)
end

function _apply_measure_instruction(instr::AbstractMeasureInstruction, state::TrackedState)
    (objs_after, results) = apply_decay_instruction(instr, state.objs)
    return (TrackedState(objs_after, state.tracker), results)
end

for T in (
    :MeasurePolar,
    :MeasureSpherical,
    :MeasureMassCosThetaPhi,
    :MeasureCosThetaPhi,
    :MeasureInvariant,
)
    @eval apply_decay_instruction(instr::$T, state::TrackedState) =
        _apply_measure_instruction(instr, state)
end

"""
    compare_instruction_paths(path1, path2, objs; T=Float64)

Run two instruction paths on the same `objs` with tracking enabled and return:
- `tracker1`, `tracker2`: accumulated Lorentz trackers for each path
- `relative`: `tracker2 * inv(tracker1)` (other relative to reference)
- `results1`, `results2`: measurement outputs along each path
- `final_objs1`, `final_objs2`: transformed objects
"""
function compare_instruction_paths(path1, path2, objs; T::Type{<:Real}=Float64)
    (state1, results1) = apply_decay_instruction(path1, init_tracked_state(objs; T = T))
    (state2, results2) = apply_decay_instruction(path2, init_tracked_state(objs; T = T))
    rel = relative_tracker(state1.tracker, state2.tracker)

    return (
        tracker1 = state1.tracker,
        tracker2 = state2.tracker,
        relative = rel,
        results1 = results1,
        results2 = results2,
        final_objs1 = state1.objs,
        final_objs2 = state2.objs,
    )
end
