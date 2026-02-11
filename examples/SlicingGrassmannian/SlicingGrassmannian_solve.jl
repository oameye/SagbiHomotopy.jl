using SagbiHomotopy
using HomotopyContinuation

# Warm-up direct solve
let
    gr = grassmannian_parametrization(2, 4)
    A = randn(ComplexF64, 4, length(gr.parametrization.groups[1]))
    S = System(A * gr.parametrization.groups[1])
    HomotopyContinuation.solve(S; show_progress = false)
end

j = 7

for m in 4:j
    for k in 2:floor(Int, m / 2)
        gr = grassmannian_parametrization(k, m)
        sagbi = gr.parametrization.groups[1]

        A = randn(ComplexF64, k * (m - k), length(sagbi))
        S = System(A * sagbi)

        result, elapsed_time, allocations, _, _ = @timed HomotopyContinuation.solve(S; show_progress = false)

        println("k=$k m=$m direct_solutions=$(length(solutions(result))) time=$(elapsed_time)s alloc=$(allocations)")
    end
end
