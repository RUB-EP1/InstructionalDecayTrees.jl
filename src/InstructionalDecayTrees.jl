module InstructionalDecayTrees

# Dependencies
using FourVectors
using LinearAlgebra

# instructions.jl
export ToHelicityFrame,
    ToHelicityFrameParticle2,
    PlaneAlign,
    ToGottfriedJacksonFrame,
    MeasurePolar,
    MeasureSpherical,
    MeasureMassCosThetaPhi,
    MeasureCosThetaPhi,
    MeasureInvariant
include("instructions.jl")

# composite_instructions.jl
export CompositeInstruction, Fork
include("composite_instructions.jl")

# execution.jl
export apply_decay_instruction,
    execute_decay_program  # Deprecated, use apply_decay_instruction instead
include("execution.jl")

# backend_fourvectors.jl
include("backend_fourvectors.jl")

# lorentz_math.jl
include("lorentz_math.jl")

# lorentz_tracker.jl
export LorentzTracker,
    relative_tracker,
    decode_lorentz_helicity,
    wigner_zyz
include("lorentz_tracker.jl")

# tracked_state.jl
export TrackedState,
    init_tracked_state,
    compare_instruction_paths
include("tracked_state.jl")

end # module InstructionalDecayTrees
