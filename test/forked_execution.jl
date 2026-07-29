struct CountStep{C} <: InstructionalDecayTrees.AbstractInstruction
    counter::C
end

function InstructionalDecayTrees.apply_decay_instruction(instr::CountStep, state)
    instr.counter[] += 1
    return (state + 1, (;))
end

struct RecordState{Tag} <: InstructionalDecayTrees.AbstractMeasureInstruction end

function InstructionalDecayTrees.apply_decay_instruction(::RecordState{Tag}, state) where {Tag}
    return (state, NamedTuple{(Tag,)}((state,)))
end

struct RecordTracker{Tag} <: InstructionalDecayTrees.AbstractMeasureInstruction end

function InstructionalDecayTrees.apply_decay_instruction(
    ::RecordTracker{Tag},
    state::TrackedState,
) where {Tag}
    return (state, NamedTuple{(Tag,)}((state.tracker,)))
end

function fourvectors_approx(xs, ys; atol=1e-12)
    return all(
        isapprox(x.E, y.E; atol = atol) &&
        isapprox(x.px, y.px; atol = atol) &&
        isapprox(x.py, y.py; atol = atol) &&
        isapprox(x.pz, y.pz; atol = atol) for (x, y) in zip(xs, ys)
    )
end

const FORK_TEST_OBJECTS = (
    FourVector(0.40, 0.10, 0.20; M = 0.30),
    FourVector(-0.20, 0.25, -0.10; M = 0.40),
    FourVector(0.05, -0.30, 0.35; M = 0.20),
    FourVector(-0.10, -0.10, -0.25; M = 0.50),
)

@testset "Fork representation and empty fork" begin
    fork = Fork((
        (MeasurePolar(:left, 1),),
        (MeasurePolar(:right, 2),),
    ))

    @test fork isa InstructionalDecayTrees.AbstractInstruction
    @test isconcretetype(typeof(fork))
    @test repr(fork) ==
          "Fork(((MeasurePolar(:left, 1),), (MeasurePolar(:right, 2),)))"

    state = init_tracked_state(FORK_TEST_OBJECTS)
    returned, results = apply_decay_instruction(Fork(()), state)
    @test returned === state
    @test isempty(results)
end

@testset "Balanced fork matches independent paths" begin
    root = ToHelicityFrame((1, 2, 3, 4))
    root_measurement = MeasureInvariant(:root_m2, (1, 2, 3, 4))
    left_branch = (
        ToHelicityFrame((1, 2)),
        MeasureCosThetaPhi(:left, 1),
    )
    right_branch = (
        ToHelicityFrameParticle2((3, 4)),
        MeasureCosThetaPhi(:right, 3),
    )
    program = (
        root,
        root_measurement,
        Fork((left_branch, right_branch)),
    )

    final_state, results = apply_decay_instruction(program, FORK_TEST_OBJECTS)
    parent_state, parent_results =
        apply_decay_instruction((root, root_measurement), FORK_TEST_OBJECTS)
    _, left_results =
        apply_decay_instruction((root, root_measurement, left_branch...), FORK_TEST_OBJECTS)
    _, right_results =
        apply_decay_instruction((root, root_measurement, right_branch...), FORK_TEST_OBJECTS)

    @test fourvectors_approx(final_state, parent_state)
    @test results.root_m2 ≈ parent_results.root_m2
    @test results.left == left_results.left
    @test results.right == right_results.right
end

@testset "Nested sequential fork returns to each parent" begin
    root = ToHelicityFrame((1, 2, 3, 4))
    to_123 = ToHelicityFrame((1, 2, 3))
    to_12 = ToHelicityFrame((1, 2))

    program = (
        root,
        Fork((
            (
                to_123,
                MeasureInvariant(:m123, (1, 2, 3)),
                Fork((
                    (
                        to_12,
                        MeasureCosThetaPhi(:p1_in_12, 1),
                    ),
                    (MeasureCosThetaPhi(:p3_in_123, 3),),
                )),
                MeasureCosThetaPhi(:after_nested, 3),
            ),
            (MeasureCosThetaPhi(:p4_at_root, 4),),
        )),
    )

    final_state, results = apply_decay_instruction(program, FORK_TEST_OBJECTS)
    parent_state, _ = apply_decay_instruction((root,), FORK_TEST_OBJECTS)
    _, path_123 = apply_decay_instruction(
        (root, to_123, MeasureInvariant(:m123, (1, 2, 3))),
        FORK_TEST_OBJECTS,
    )
    _, path_12 = apply_decay_instruction(
        (root, to_123, to_12, MeasureCosThetaPhi(:p1_in_12, 1)),
        FORK_TEST_OBJECTS,
    )
    _, path_4 = apply_decay_instruction(
        (root, MeasureCosThetaPhi(:p4_at_root, 4)),
        FORK_TEST_OBJECTS,
    )

    @test fourvectors_approx(final_state, parent_state)
    @test results.m123 ≈ path_123.m123
    @test results.p1_in_12 == path_12.p1_in_12
    @test results.p4_at_root == path_4.p4_at_root
    @test results.after_nested == results.p3_in_123
end

@testset "Tracked branches inherit and restore the parent tracker" begin
    root = ToHelicityFrame((1, 2, 3, 4))
    left = ToHelicityFrame((1, 2))
    right = ToHelicityFrameParticle2((3, 4))
    initial = init_tracked_state(FORK_TEST_OBJECTS)

    program = (
        root,
        Fork((
            (left, RecordTracker{:left_tracker}()),
            (right, RecordTracker{:right_tracker}()),
        )),
    )
    final_state, results = apply_decay_instruction(program, initial)
    parent_state, _ = apply_decay_instruction((root,), initial)
    left_state, _ = apply_decay_instruction((root, left), initial)
    right_state, _ = apply_decay_instruction((root, right), initial)

    @test final_state.tracker.Λ ≈ parent_state.tracker.Λ atol = 1e-12
    @test final_state.tracker.U ≈ parent_state.tracker.U atol = 1e-12
    @test results.left_tracker.Λ ≈ left_state.tracker.Λ atol = 1e-12
    @test results.left_tracker.U ≈ left_state.tracker.U atol = 1e-12
    @test results.right_tracker.Λ ≈ right_state.tracker.Λ atol = 1e-12
    @test results.right_tracker.U ≈ right_state.tracker.U atol = 1e-12
end

@testset "Fork rejects duplicate measurement tags" begin
    @test_throws ArgumentError apply_decay_instruction(
        Fork(((RecordState{:x}(),), (RecordState{:x}(),))),
        0,
    )
    @test_throws ArgumentError apply_decay_instruction(
        (RecordState{:x}(), Fork(((RecordState{:x}(),),))),
        0,
    )
    @test_throws ArgumentError apply_decay_instruction(
        (Fork(((RecordState{:x}(),),)), RecordState{:x}()),
        0,
    )
    @test_throws ArgumentError apply_decay_instruction(
        Fork(((RecordState{:x}(), RecordState{:x}()),)),
        0,
    )
    @test_throws ArgumentError apply_decay_instruction(
        (
            RecordState{:x}(),
            CompositeInstruction((Fork(((RecordState{:x}(),),)),)),
        ),
        0,
    )

    # Preserve the historical merge behavior for programs that are entirely linear.
    _, linear_results =
        apply_decay_instruction((RecordState{:x}(), RecordState{:x}()), 0)
    @test linear_results == (x = 0,)
end

@testset "Static fork inference, allocation, and shared-prefix execution" begin
    counter = Ref(0)
    program = (
        CountStep(counter),
        Fork((
            (RecordState{:left}(),),
            (RecordState{:right}(),),
        )),
    )

    branch_state, branch_results = @inferred apply_decay_instruction(program[2], 1)
    @test branch_state == 1
    @test branch_results == (left = 1, right = 1)

    final_state, results = apply_decay_instruction(program, 0)
    @test final_state == 1
    @test results == (left = 1, right = 1)
    @test counter[] == 1

    apply_decay_instruction(program, 0) # warm up the exact call before measuring
    counter[] = 0
    allocated = @allocated apply_decay_instruction(program, 0)
    @test counter[] == 1
    @test allocated <= 512
end
