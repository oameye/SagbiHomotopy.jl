using BenchmarkTools
using Random
using SagbiHomotopy
using HomotopyContinuation
using MultivariatePolynomials: subs

function oscillator_instance(N::Int, M::Int; seed::Int = 20260210)
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

    w = vcat(fill(2, M * N), fill(1, M * N))
    full_eqs = [subs(lin_sys[i][1], collect(z[i, :]) => sagbi[i]) for i in 1:length(sagbi)]

    return (
        linear = LinearSection(lin_sys),
        parametrization = Parametrization(sagbi),
        w = w,
        full_eqs = full_eqs,
    )
end

function benchmark_oscillators!(suite::BenchmarkTools.BenchmarkGroup)
    grp = BenchmarkTools.BenchmarkGroup()
    suite["Oscillators"] = grp

    inst_12 = oscillator_instance(1, 2; seed = 20260210)
    case_12 = BenchmarkTools.BenchmarkGroup()
    grp["N=1 M=2"] = case_12
    case_12["solve"] = @benchmarkable SagbiHomotopy.solve(
        $(inst_12.linear),
        $(inst_12.parametrization);
        weight = $(inst_12.w),
        check_degree = false,
    ) seconds = 20
    case_12["direct solve"] = @benchmarkable HomotopyContinuation.solve(
        $(inst_12.full_eqs);
        show_progress = false,
    ) seconds = 20

    inst_22 = oscillator_instance(2, 2; seed = 20260211)
    case_22 = BenchmarkTools.BenchmarkGroup()
    grp["N=2 M=2"] = case_22
    case_22["solve"] = @benchmarkable SagbiHomotopy.solve(
        $(inst_22.linear),
        $(inst_22.parametrization);
        weight = $(inst_22.w),
        check_degree = false,
    ) seconds = 20
    case_22["direct solve"] = @benchmarkable HomotopyContinuation.solve(
        $(inst_22.full_eqs);
        show_progress = false,
    ) seconds = 20

    return nothing
end
