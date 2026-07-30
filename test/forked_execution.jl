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

@testset "Fork construction and execution" begin
    root = ToHelicityFrame((1, 2, 3, 4))
    left = (
        ToHelicityFrame((1, 2)),
        MeasureCosThetaPhi(:left, 1),
    )
    right = (
        ToHelicityFrameParticle2((3, 4)),
        MeasureCosThetaPhi(:right, 3),
    )
    fork = Fork((left, right))

    @test fork isa InstructionalDecayTrees.AbstractInstruction
    @test isconcretetype(typeof(fork))
    @test repr(Fork(())) == "Fork(())"

    final_state, results =
        apply_decay_instruction((root, fork), FORK_TEST_OBJECTS)
    parent_state, _ =
        apply_decay_instruction((root,), FORK_TEST_OBJECTS)
    _, left_results =
        apply_decay_instruction((root, left...), FORK_TEST_OBJECTS)
    _, right_results =
        apply_decay_instruction((root, right...), FORK_TEST_OBJECTS)

    @test fourvectors_approx(final_state, parent_state)
    @test results.left == left_results.left
    @test results.right == right_results.right

    returned, empty_results =
        apply_decay_instruction(Fork(()), parent_state)
    @test returned === parent_state
    @test isempty(empty_results)
end

@testset "Nested and tracked forks return to their parent" begin
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

    final_state, results =
        apply_decay_instruction(program, FORK_TEST_OBJECTS)
    parent_state, _ =
        apply_decay_instruction((root,), FORK_TEST_OBJECTS)
    _, path_12 = apply_decay_instruction(
        (root, to_123, to_12, MeasureCosThetaPhi(:p1_in_12, 1)),
        FORK_TEST_OBJECTS,
    )
    _, path_4 = apply_decay_instruction(
        (root, MeasureCosThetaPhi(:p4_at_root, 4)),
        FORK_TEST_OBJECTS,
    )

    @test fourvectors_approx(final_state, parent_state)
    @test results.p1_in_12 == path_12.p1_in_12
    @test results.p4_at_root == path_4.p4_at_root
    @test results.after_nested == results.p3_in_123

    tracked = init_tracked_state(FORK_TEST_OBJECTS)
    final_tracked, tracked_results =
        apply_decay_instruction(program, tracked)
    parent_tracked, _ =
        apply_decay_instruction((root,), tracked)

    @test final_tracked.tracker.Λ ≈ parent_tracked.tracker.Λ atol = 1e-12
    @test final_tracked.tracker.U ≈ parent_tracked.tracker.U atol = 1e-12
    @test tracked_results == results
end

@testset "Fork rejects duplicate measurement tags" begin
    left = (MeasurePolar(:angle, 1),)
    right = (MeasurePolar(:angle, 2),)

    @test_throws ArgumentError apply_decay_instruction(
        Fork((left, right)),
        FORK_TEST_OBJECTS,
    )
    @test_throws ArgumentError apply_decay_instruction(
        (MeasurePolar(:angle, 1), Fork((right,))),
        FORK_TEST_OBJECTS,
    )
    @test_throws ArgumentError apply_decay_instruction(
        Fork(((MeasurePolar(:angle, 1), MeasurePolar(:angle, 2)),)),
        FORK_TEST_OBJECTS,
    )

    # Linear programs retain their historical later-value-wins behavior.
    _, linear_results = apply_decay_instruction(
        (MeasurePolar(:angle, 1), MeasurePolar(:angle, 2)),
        FORK_TEST_OBJECTS,
    )
    _, expected =
        apply_decay_instruction(MeasurePolar(:angle, 2), FORK_TEST_OBJECTS)
    @test linear_results == expected
end

@testset "Static empty-result fork inference" begin
    fork = Fork((
        (ToHelicityFrame((1, 2, 3, 4)),),
        (),
    ))

    state, results =
        @inferred apply_decay_instruction(fork, FORK_TEST_OBJECTS)
    @test state === FORK_TEST_OBJECTS
    @test isempty(results)

    apply_decay_instruction(fork, FORK_TEST_OBJECTS)
    @test (@allocated apply_decay_instruction(fork, FORK_TEST_OBJECTS)) <= 512
end
