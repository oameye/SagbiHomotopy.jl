using SagbiHomotopy
using HomotopyContinuation
using Random
using Test

function _nonzero_int_matrix(rng::AbstractRNG, m::Int, n::Int)
    values = collect(-100:-1)
    append!(values, 1:100)
    return rand(rng, values, m, n)
end

@testset "V2 API shape" begin
    @var x y z

    p1 = Parametrization([x, y, z])
    @test length(p1.groups) == 1
    @test p1.groups[1] == [x, y, z]

    p2 = Parametrization([[x, y], [y, z]])
    @test length(p2.groups) == 2

    l1 = LinearSection([x + y, y + z])
    @test length(l1.groups) == 1

    l2 = LinearSection([[x + y, y + z], [x + z]])
    @test length(l2.groups) == 2

    @test_throws InputShapeError Parametrization(Any[[[x]], [[y]]])
    @test_throws InputShapeError LinearSection(Any[[[x + y]]])

    @var p[1:2]
    @var q[1:3]
    linear = LinearSection([p[1] + p[2]])
    param = Parametrization([q[1], q[2], q[3]])
    @test_throws IncompatibleProblemError SagbiProblem(linear, param)
end

@testset "V2 errors and determinism" begin
    @var x y
    @test_throws NoSagbiWeightError detect_weight(Parametrization([x + y, x - y]))

    @var a[1:4]
    sagbi = [[x * (x^2 + y^2 - 2 * x), x * (5 - 4 * y), y * (x^2 + y^2 - 2 * x), y * (5 - 4 * y)]]
    weight = Int[-2, -1]
    linear = LinearSection([a[1] + 2a[2] + 3a[3] + 4a[4], 2a[1] + 5a[2] + 7a[3] + 11a[4]])
    param = Parametrization(sagbi)

    @test_throws DegreeDropError SagbiHomotopy.solve(linear, param; weight = weight, check_degree = true)

    @var u v w
    param_det = Parametrization([u^2 + 1, v^2 + 1, u * v + w^2, 1])
    weight_det = Int[-2, -2, -3]

    rng1 = MersenneTwister(12345)
    rng2 = MersenneTwister(12345)
    @test degree_map(param_det; rng = rng1) == degree_map(param_det; rng = rng2)

    @var r[1:4]
    A = _nonzero_int_matrix(MersenneTwister(2026), 3, 4)
    linear_det = LinearSection(A * r)

    opts1 = SolveOptions(weight = weight_det, check_degree = true, allow_degree_drop = true, rng = MersenneTwister(17))
    opts2 = SolveOptions(weight = weight_det, check_degree = true, allow_degree_drop = true, rng = MersenneTwister(17))
    res1 = SagbiHomotopy.solve(linear_det, param_det, opts1)
    res2 = SagbiHomotopy.solve(linear_det, param_det, opts2)

    @test res1.metadata.degree_map == res2.metadata.degree_map
    @test res1.metadata.degree_monomial_map == res2.metadata.degree_monomial_map
end

@testset "V2 inferred constructors" begin
    @var x y z
    @test @inferred(Parametrization([x, y, z])) isa Parametrization
    @test @inferred(LinearSection([x + y, y + z])) isa LinearSection

    p = @inferred Parametrization([x^2 + 1, y^2 + 1, x * y + z^2, 1])
    @test @inferred(degree_map(p; rng = MersenneTwister(1))) == 8
end
