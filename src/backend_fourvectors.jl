# --- Helper Functions ---

"""
    get_fourvector(objs, indices::Tuple)

Extract a four-vector from `objs` based on `indices`:
- Single positive index: returns `objs[index]`
- Single negative index: returns `-objs[-index]`
- Multiple indices: returns sum with each index handled independently
  (e.g., `(1, -2, 3)` → `objs[1] + (-objs[2]) + objs[3]`)
"""
function get_fourvector(objs, indices::Tuple{Vararg{Int}})
    if length(indices) == 1
        idx = first(indices)
        return idx < 0 ? -objs[-idx] : objs[idx]
    else
        return sum(idx < 0 ? -objs[-idx] : objs[idx] for idx in indices)
    end
end

# --- Implementations ---

function apply_decay_instruction(instr::ToHelicityFrame, objs)
    # Get four-vector from indices (handles single, multiple, and negative indices)
    P_tot = get_fourvector(objs, instr.indices)

    # Apply boost to all
    new_objs = map(p -> transform_to_cmf(p, P_tot), objs)

    return (new_objs, _empty_instruction_results)
end

function apply_decay_instruction(instr::ToHelicityFrameParticle2, objs)
    P_tot = get_fourvector(objs, instr.indices)

    # Uses -vec for Rz and Ry
    # Construct vector with opposite spatial components to get angles
    P_inv = FourVector(-P_tot.px, -P_tot.py, -P_tot.pz; E = P_tot.E)

    ϕ_inv = azimuthal_angle(P_inv)
    θ_inv = polar_angle(P_inv)
    γ = boost_gamma(P_tot)

    # Sequence: Rz(-ϕ_inv) -> Ry(-θ_inv) -> Ry(-π) -> Bz(-γ)
    transform = p -> p |> Rz(-ϕ_inv) |> Ry(-θ_inv) |> Ry(-π) |> Bz(-γ)

    new_objs = map(transform, objs)
    return (new_objs, _empty_instruction_results)
end

function apply_decay_instruction(instr::PlaneAlign, objs)
    # Get four-vectors from index specifications (handles single, multiple, and negative indices)
    axis_z = get_fourvector(objs, instr.z_idx)
    axis_x = get_fourvector(objs, instr.x_idx)

    final_objs = map(p -> rotate_to_plane(p, axis_z, axis_x), objs)

    return (final_objs, _empty_instruction_results)
end

# CompositeInstruction handling is now in execution.jl
# This method is no longer needed since apply_decay_instruction is defined there

function apply_decay_instruction(instr::ToGottfriedJacksonFrame, objs)
    # Step 1: Apply ToHelicityFrame (reuses existing implementation)
    hel_instr = ToHelicityFrame(instr.system_indices)
    (objs_after_boost, _) = apply_decay_instruction(hel_instr, objs)

    # Step 2: Apply PlaneAlign (reuses existing implementation)
    # Standard GJ definition: target_idx is always positive but negated when passed to PlaneAlign
    # Negate each index in target_idx tuple
    target_idx_negated = Tuple(-idx for idx in instr.target_idx)
    plane_instr = PlaneAlign(instr.beam_idx, target_idx_negated)
    (final_objs, _) = apply_decay_instruction(plane_instr, objs_after_boost)

    return (final_objs, _empty_instruction_results)
end

function apply_decay_instruction(instr::MeasurePolar, objs)
    p = get_fourvector(objs, instr.idx)
    val = polar_angle(p)
    return (objs, NamedTuple{(instr.tag,)}((val,)))
end

function apply_decay_instruction(instr::MeasureSpherical, objs)
    p = get_fourvector(objs, instr.indices)
    val_theta = polar_angle(p)
    val_phi = azimuthal_angle(p)
    # Construct result tuple with both keys
    res = NamedTuple{(instr.theta_tag, instr.phi_tag)}((val_theta, val_phi))
    return (objs, res)
end

function apply_decay_instruction(instr::MeasureMassCosThetaPhi, objs)
    p = get_fourvector(objs, instr.indices)

    m_val = mass(p)
    cos_theta_val = cos_theta(p)
    phi_val = azimuthal_angle(p)

    # Structure: tag => (m, cosθ, ϕ)
    val_tuple = (m = m_val, cosθ = cos_theta_val, ϕ = phi_val)

    return (objs, NamedTuple{(instr.tag,)}((val_tuple,)))
end

function apply_decay_instruction(instr::MeasureCosThetaPhi, objs)
    p = get_fourvector(objs, instr.indices)

    cos_theta_val = cos_theta(p)
    phi_val = azimuthal_angle(p)

    # Structure: tag => (cosθ, ϕ)
    val_tuple = (cosθ = cos_theta_val, ϕ = phi_val)

    return (objs, NamedTuple{(instr.tag,)}((val_tuple,)))
end

function apply_decay_instruction(instr::MeasureInvariant, objs)
    P_tot = get_fourvector(objs, instr.indices)
    val = mass(P_tot)^2
    return (objs, NamedTuple{(instr.tag,)}((val,)))
end
