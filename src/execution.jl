"""
    apply_decay_instruction(instr, objs)

Execute an instruction or instruction sequence on `objs`.
Returns a tuple: `(modified_objects, results_named_tuple)`.

This is the public execution entry point. The `instr` can be:
- A single instruction: Executed directly
- A `CompositeInstruction`: Execute a reusable instruction group
- A tuple of instructions: treated as a convenient instruction sequence
- A `Fork`: execute branches from one parent state and return to that state

Every method returns `(state_or_objs, measurement_results)`. Frame transforms
produce no measurements, so they use `_empty_instruction_results`.
"""
function apply_decay_instruction end

# Empty measurement results returned by frame-transform instructions.
const _empty_instruction_results = (;)

function _merge_instruction_results(left::NamedTuple, right::NamedTuple, ::Val{false})
    return merge(left, right)
end

_apply_instruction(instr, state, ::Val) = apply_decay_instruction(instr, state)

function _apply_instruction_sequence(
    instructions::Tuple,
    state,
    mode::Val,
)
    current_state = state
    all_results = _empty_instruction_results
    for instruction in instructions
        current_state, instruction_results =
            _apply_instruction(instruction, current_state, mode)
        all_results = _merge_instruction_results(all_results, instruction_results, mode)
    end
    return (current_state, all_results)
end

function apply_decay_instruction(instructions::Tuple, state)
    return apply_decay_instruction(CompositeInstruction(instructions), state)
end

function apply_decay_instruction(instr::CompositeInstruction, state)
    mode = Val(_contains_fork(instr))
    return _apply_instruction_sequence(instr.instructions, state, mode)
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
