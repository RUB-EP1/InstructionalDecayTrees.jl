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

"""
    fork_branch_state(state)

Create the independent state snapshot used to start one [`Fork`](@ref) branch.

The default uses `deepcopy` so mutable state in one branch cannot affect its
siblings or the parent returned after the fork. Backends whose states are
immutable or persistent may overload this function and safely return `state`
itself to avoid the copy.
"""
fork_branch_state(state) = deepcopy(state)

_contains_fork(::Any) = false
_contains_fork(::Fork) = true
_contains_fork(instr::CompositeInstruction) = _contains_fork(instr.instructions)

function _contains_fork(instructions::Tuple)
    for instruction in instructions
        _contains_fork(instruction) && return true
    end
    return false
end

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
_apply_instruction(instr::Tuple, state, mode::Val{true}) =
    _apply_instruction_sequence(instr, state, mode)
_apply_instruction(instr::CompositeInstruction, state, mode::Val{true}) =
    _apply_instruction_sequence(instr.instructions, state, mode)
_apply_instruction(instr::Fork, state, ::Val) = _apply_fork(instr, state)

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
    branch_state = fork_branch_state(parent_state)
    _, branch_results = _apply_instruction(first(branches), branch_state, Val(true))
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
