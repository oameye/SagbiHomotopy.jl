using SagbiHomotopy
using HomotopyContinuation

# Warm-up
let
    gr = grassmannian_parametrization(2, 4)
    @var z[1:length(gr.parametrization.groups[1])]
    A = randn(ComplexF64, 4, length(z))
    linear = LinearSection(A * z)
    SagbiHomotopy.solve(linear, gr.parametrization; weight = gr.weight, check_degree = false)
end

j = 7

for m in 4:j
    for k in 2:floor(Int, m / 2)
        gr = grassmannian_parametrization(k, m)
        ncoords = length(gr.parametrization.groups[1])

        @var z[1:ncoords]
        A = randn(ComplexF64, k * (m - k), ncoords)
        linear = LinearSection(A * z)

        result, elapsed_time, allocations, _, _ = @timed SagbiHomotopy.solve(
            linear,
            gr.parametrization;
            weight = gr.weight,
            check_degree = false,
        )

        println("k=$k m=$m solutions=$(length(result.solutions)) time=$(elapsed_time)s alloc=$(allocations)")
    end
end
