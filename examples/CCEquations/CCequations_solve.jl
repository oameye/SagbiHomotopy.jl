using SagbiHomotopy
using HomotopyContinuation

# Warm-up direct solve
let
    gr = grassmannian_parametrization(2, 4)
    @var l
    A = randn(ComplexF64, 5, length(gr.parametrization.groups[1]))
    equations = A * gr.parametrization.groups[1] - l .* vcat(1, gr.vars)
    HomotopyContinuation.solve(equations; show_progress = false)
end

j = 8

for m in 4:j
    for k in 2:floor(Int, m / 2)
        gr = grassmannian_parametrization(k, m)
        @var l

        A = randn(ComplexF64, k * (m - k) + 1, length(gr.parametrization.groups[1]))
        equations = A * gr.parametrization.groups[1] - l .* vcat(1, gr.vars)

        result, elapsed_time, allocations, _, _ = @timed HomotopyContinuation.solve(equations; show_progress = false)

        println("k=$k m=$m direct_solutions=$(length(solutions(result))) time=$(elapsed_time)s alloc=$(allocations)")
    end
end
