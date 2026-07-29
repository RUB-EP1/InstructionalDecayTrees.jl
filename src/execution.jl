"""
    apply_decay_instruction(instr, objs)

Execute an instruction or instruction sequence on `objs`.
Returns a tuple: `(modified_objects, results_named_tuple)`.

This is the public execution entry point. The `instr` can be:
- A single instruction: Executed directly
- A `CompositeInstruction`: Executed with nested recursive execution, maintaining encapsulation
- A tuple of instructions: treated as a convenient instruction sequence
- A `Fork`: execute branches from one parent state and return to that state

Nested `CompositeInstruction`s are executed recursively, maintaining encapsulation
of reusable instruction groups.

Every method returns `(state_or_objs, measurement_results)`. Frame transforms
produce no measurements, so they use `_empty_instruction_results`.
"""
function apply_decay_instruction end

# Empty measurement results returned by frame-transform instructions.
const _empty_instruction_results = (;)

# A fork makes result merging strict for its entire enclosing instruction tree.
# Carry the policy as a `Val` so statically shaped tuples remain specialized.
@inline _fork_mode(::Any) = Val(false)
@inline _fork_mode(::Fork) = Val(true)
@inline _fork_mode(instr::CompositeInstruction) = _fork_mode(instr.instructions)
@inline _fork_mode(::Tuple{}) = Val(false)
@inline _fork_mode(instrs::Tuple) =
    _combine_fork_modes(_fork_mode(first(instrs)), _fork_mode(Base.tail(instrs)))

@inline _combine_fork_modes(::Val{true}, ::Val) = Val(true)
@inline _combine_fork_modes(::Val{false}, tail_mode::Val) = tail_mode

function _merge_instruction_results(left::NamedTuple, right::NamedTuple, ::Val{false})
    return merge(left, right)
end

function _merge_instruction_results(left::NamedTuple, right::NamedTuple, ::Val{true})
    for tag in keys(right)
        haskey(left, tag) && throw(ArgumentError(
            "duplicate measurement tag $(repr(tag)) in a forked instruction program",
        ))
    end
    return merge(left, right)
end

_apply_instruction(instr, state, ::Val) = apply_decay_instruction(instr, state)
_apply_instruction(instr::Tuple, state, mode::Val) =
    _apply_instruction_sequence(instr, state, mode)
_apply_instruction(instr::CompositeInstruction, state, mode::Val) =
    _apply_instruction_sequence(instr.instructions, state, mode)
_apply_instruction(instr::Fork, state, ::Val) = _apply_fork(instr, state)

function _apply_instruction_sequence(
    ::Tuple{},
    state,
    ::Val,
    results::NamedTuple=_empty_instruction_results,
)
    return (state, results)
end

function _apply_instruction_sequence(
    instructions::Tuple,
    state,
    mode::Val,
    results::NamedTuple=_empty_instruction_results,
)
    instruction = first(instructions)
    next_state, instruction_results = _apply_instruction(instruction, state, mode)
    next_results = _merge_instruction_results(results, instruction_results, mode)
    return _apply_instruction_sequence(
        Base.tail(instructions),
        next_state,
        mode,
        next_results,
    )
end

function apply_decay_instruction(instructions::Tuple, state)
    mode = _fork_mode(instructions)
    return _apply_instruction_sequence(instructions, state, mode)
end

function apply_decay_instruction(instr::CompositeInstruction, state)
    mode = _fork_mode(instr)
    return _apply_instruction_sequence(instr.instructions, state, mode)
end

function _apply_fork_branches(
    ::Tuple{},
    parent_state,
    results::NamedTuple=_empty_instruction_results,
)
    return results
end

function _apply_fork_branches(
    branches::Tuple,
    parent_state,
    results::NamedTuple=_empty_instruction_results,
)
    _, branch_results = _apply_instruction(first(branches), parent_state, Val(true))
    next_results = _merge_instruction_results(results, branch_results, Val(true))
    return _apply_fork_branches(Base.tail(branches), parent_state, next_results)
end

function _apply_fork(instr::Fork, state)
    results = _apply_fork_branches(instr.branches, state)
    return (state, results)
end

apply_decay_instruction(instr::Fork, state) = _apply_fork(instr, state)

"""
    execute_decay_program(objs, program)

Deprecated: use `apply_decay_instruction(sequence, objs)` instead.
This function is kept for backward compatibility but will be removed in a future version.
"""
function execute_decay_program(objs, program)
    Base.depwarn("`execute_decay_program` is deprecated, use `apply_decay_instruction` instead", :execute_decay_program)
    return apply_decay_instruction(program, objs)
end
