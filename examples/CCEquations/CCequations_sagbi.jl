using SagbiHomotopy
using HomotopyContinuation
using MultivariatePolynomials: subs

# Warm-up
let
    gr = grassmannian_parametrization(2, 4)
    @var l
    basis = vcat(l .* vcat(1, gr.vars), gr.parametrization.groups[1])
    @var z[1:length(basis)]

    A = randn(ComplexF64, 5, length(gr.parametrization.groups[1]))
    equations = A * gr.parametrization.groups[1] - l .* vcat(1, gr.vars)
    linear = LinearSection(subs(equations, basis => z))

    param = Parametrization(basis)
    w = vcat(1, gr.weight)
    SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
end

j = 8

for m in 4:j
    for k in 2:floor(Int, m / 2)
        gr = grassmannian_parametrization(k, m)
        @var l

        basis = vcat(l .* vcat(1, gr.vars), gr.parametrization.groups[1])
        @var z[1:length(basis)]

        A = randn(ComplexF64, k * (m - k) + 1, length(gr.parametrization.groups[1]))
        equations = A * gr.parametrization.groups[1] - l .* vcat(1, gr.vars)

        linear = LinearSection(subs(equations, basis => z))
        param = Parametrization(basis)
        w = vcat(1, gr.weight)

        result, elapsed_time, allocations, _, _ = @timed SagbiHomotopy.solve(
            linear,
            param;
            weight = w,
            check_degree = false,
        )

        println("k=$k m=$m solutions=$(length(result.solutions)) time=$(elapsed_time)s alloc=$(allocations)")
    end
end
