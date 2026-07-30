"""
    fork_branch_state(state)

Create the independent state snapshot used to start one [`Fork`](@ref) branch.

This is an extension point for execution backends; application code normally
does not call it directly. The default uses `deepcopy` so mutable state in one
branch cannot affect its siblings or the parent returned after the fork.
Backends with immutable or persistent states may overload this method to return
`state` itself.
"""
fork_branch_state(state) = deepcopy(state)

_contains_fork(::Any) = false
_contains_fork(::Fork) = true
_contains_fork(instr::CompositeInstruction) = _contains_fork(instr.instructions)

function _contains_fork(instructions::Tuple)
    return any(_contains_fork, instructions)
end

function _merge_instruction_results(left::NamedTuple, right::NamedTuple, ::Val{true})
    for tag in keys(right)
        haskey(left, tag) && throw(ArgumentError(
            "duplicate measurement tag $(repr(tag)) in a forked instruction program",
        ))
    end
    return merge(left, right)
end

# A program containing a fork uses strict result merging throughout its nested
# tuples and composites, so duplicate tags cannot hide at a nesting boundary.
_apply_instruction(instr::Tuple, state, mode::Val{true}) =
    _apply_instruction_sequence(instr, state, mode)
_apply_instruction(instr::CompositeInstruction, state, mode::Val{true}) =
    _apply_instruction_sequence(instr.instructions, state, mode)

_fork_results(::Tuple{}, parent_state) = _empty_instruction_results

function _fork_results(branches::Tuple, parent_state)
    branch_state = fork_branch_state(parent_state)
    _, branch_results =
        _apply_instruction(first(branches), branch_state, Val(true))
    remaining_results = _fork_results(Base.tail(branches), parent_state)
    return _merge_instruction_results(branch_results, remaining_results, Val(true))
end

function apply_decay_instruction(instr::Fork, state)
    return (state, _fork_results(instr.branches, state))
end
