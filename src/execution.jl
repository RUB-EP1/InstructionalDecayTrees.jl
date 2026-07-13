"""
    apply_decay_instruction(instr, objs)

Execute an instruction or instruction sequence on `objs`.
Returns a tuple: `(modified_objects, results_named_tuple)`.

This is the public execution entry point. The `instr` can be:
- A single instruction: Executed directly
- A `CompositeInstruction`: Executed with nested recursive execution, maintaining encapsulation
- A tuple of instructions: treated as a convenient instruction sequence

Nested `CompositeInstruction`s are executed recursively, maintaining encapsulation
of reusable instruction groups.

Every method returns `(state_or_objs, measurement_results)`. Frame transforms
produce no measurements, so they use `_empty_instruction_results`.
"""
function apply_decay_instruction end

# Empty measurement results returned by frame-transform instructions.
const _empty_instruction_results = (;)

# Convenience: accept tuples and wrap them in CompositeInstruction
function apply_decay_instruction(instr::Tuple, objs)
    return apply_decay_instruction(CompositeInstruction(instr), objs)
end

# Base case: empty composite instruction
function apply_decay_instruction(instr::CompositeInstruction{<:Tuple{}}, objs)
    return (objs, _empty_instruction_results)
end

# Execute composite instruction: iterate over tuple, recurse only for nested composites
function apply_decay_instruction(instr::CompositeInstruction, objs)
    instructions = instr.instructions
    current_objs = objs
    all_results = NamedTuple()

    # Iterate over all instructions in the tuple (no recursion needed here)
    for instruction in instructions
        # Only recurse if we hit a nested CompositeInstruction
        (current_objs, instruction_results) = apply_decay_instruction(instruction, current_objs)
        all_results = merge(all_results, instruction_results)
    end

    return (current_objs, all_results)
end

"""
    execute_decay_program(objs, program)

Deprecated: use `apply_decay_instruction(sequence, objs)` instead.
This function is kept for backward compatibility but will be removed in a future version.
"""
function execute_decay_program(objs, program)
    Base.depwarn("`execute_decay_program` is deprecated, use `apply_decay_instruction` instead", :execute_decay_program)
    return apply_decay_instruction(program, objs)
end
