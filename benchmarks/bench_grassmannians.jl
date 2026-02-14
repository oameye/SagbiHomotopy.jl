using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation

function benchmark_grassmannians!(suite::BenchmarkTools.BenchmarkGroup)
    rng = MersenneTwister(20260210)
    gr = grassmannian_parametrization(2, 6)

    @var z[1:15]
    A = randn(rng, ComplexF64, 8, 15)
    linear = LinearSection(A * z)

    suite["Grassmannians"]["Gr(2,6) slice solve"] = @benchmarkable SagbiHomotopy.solve(
        $linear,
        $(gr.parametrization);
        weight = $(gr.weight),
        check_degree = false,
    ) seconds = 20

    return nothing
end
