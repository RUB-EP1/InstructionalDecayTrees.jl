"""
    CompositeInstruction(instructions)

Group a tuple of instructions into a reusable instruction sequence.

Plain tuples of instructions can be passed directly to [`apply_decay_instruction`](@ref);
`CompositeInstruction` is useful when you want to name or dispatch on a sequence.
"""
struct CompositeInstruction{T<:Tuple} <: AbstractInstruction
    instructions::T  # Tuple of AbstractInstruction objects
end

function Base.show(io::IO, instr::CompositeInstruction)
    print(io, "CompositeInstruction(")
    show(io, instr.instructions)
    print(io, ")")
end

"""
    Fork(branches)

Execute several instruction `branches` independently from the same incoming
state and merge their measurement results.

`branches` is a tuple whose elements are instructions, instruction tuples, or
[`CompositeInstruction`](@ref)s. Branches run in tuple order. Transformations
inside one branch are never visible to its siblings, and the state returned by
the `Fork` is the unchanged incoming parent state. Forks may be nested.

Measurement tags must be unique throughout any program containing a `Fork`;
duplicate tags raise an `ArgumentError`. The tuple representation is retained
in the type so the compiler can specialize on a statically shaped tree.
"""
struct Fork{B<:Tuple} <: AbstractInstruction
    branches::B
end

function Base.show(io::IO, instr::Fork)
    print(io, "Fork(")
    show(io, instr.branches)
    print(io, ")")
end
