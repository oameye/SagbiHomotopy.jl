using SagbiHomotopy
using HomotopyContinuation
using Random
using MultivariatePolynomials: subs

function oscillator_full_system(N::Int, M::Int; seed::Int = 20260210)
    rng = MersenneTwister(seed)

    @var u[1:M, 1:N] v[1:M, 1:N]
    L = M + 3 + N - 1
    @var z[1:(2 * M * N), 1:L]

    polynomials = Matrix{HomotopyContinuation.Expression}(undef, 2 * M * N, L)
    for k in 1:N
        for i in 1:M
            r1 = 2 * i - 1 + 2 * M * (k - 1)
            r2 = 2 * i + 2 * M * (k - 1)

            polynomials[r1, 1] = (u[1, 1])^0
            polynomials[r1, 2] = u[i, k]
            polynomials[r1, 3] = v[i, k]

            polynomials[r2, 1] = (u[1, 1])^0
            polynomials[r2, 2] = u[i, k]
            polynomials[r2, 3] = v[i, k]

            for j in 1:M
                polynomials[r1, j + 3] = u[i, k] * (u[j, k]^2 + v[j, k]^2)
                polynomials[r2, j + 3] = v[i, k] * (u[j, k]^2 + v[j, k]^2)
            end

            for j in 1:(N - 1)
                lst = filter(jj -> jj != k, collect(1:N))
                polynomials[r1, j + M + 3] = v[i, lst[j]]
                polynomials[r2, j + M + 3] = u[i, lst[j]]
            end
        end
    end

    sagbi = [polynomials[i, :] for i in 1:(2 * M * N)]
    c = rand(rng, Float64, 2 * M * N, L)
    lin_sys = [[sum(z[i, :] .* c[i, :])] for i in 1:(2 * M * N)]

    full_eqs = [subs(lin_sys[i][1], collect(z[i, :]) => sagbi[i]) for i in 1:length(sagbi)]
    return full_eqs
end

# Warm-up direct solve
HomotopyContinuation.solve(oscillator_full_system(2, 1); show_progress = false)

N = 2
M = 3

for s in 1:N
    for r in 1:M
        equations = oscillator_full_system(s, r; seed = 20260210 + s + 10r)
        result, elapsed_time, allocations, _, _ = @timed HomotopyContinuation.solve(equations; show_progress = false)

        println("N=$s M=$r direct_solutions=$(length(solutions(result))) time=$(elapsed_time)s alloc=$(allocations)")
    end
end
