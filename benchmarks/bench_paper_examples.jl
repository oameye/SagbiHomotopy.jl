using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation

function benchmark_paper_examples!(suite::BenchmarkTools.BenchmarkGroup)
    @var x y z
    sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]

    param = Parametrization(sagbi)
    w = detect_weight(param)

    rng = MersenneTwister(20260210)
    @var p[1:4]
    @var q[1:4]
    lin_sys = [
        rand(rng, -100:-1, 2, 4) * p,
        rand(rng, -100:-1, 1, 4) * q,
    ]
    linear = LinearSection(lin_sys)

    suite["semimixed"]["detect_weight"] = @benchmarkable detect_weight($param) seconds = 10
    suite["semimixed"]["solve"] = @benchmarkable SagbiHomotopy.solve(
        $linear,
        $param;
        weight = $w,
        check_degree = false,
    ) seconds = 10

    return nothing
end
