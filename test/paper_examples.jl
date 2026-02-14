using SagbiHomotopy
using HomotopyContinuation
using Random
using Test

function _nonzero_int_matrix_paper(rng::AbstractRNG, m::Int, n::Int)
    values = collect(-100:-1)
    append!(values, 1:100)
    return rand(rng, values, m, n)
end

@testset "Paper examples" begin
    Random.seed!(20260210)

    @testset "Weight deformation sanity" begin
        @var x y z
        sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
        param = Parametrization(sagbi)
        w = detect_weight(param)

        @test !isempty(w)

        @unique_var t
        for poly in reduce(vcat, param.groups)
            deformed = weight_deformation(poly, param, w; t = t)
            @test subs(deformed, t => 1) == poly
            @test subs(deformed, t => 0) == SagbiHomotopy.leading_monomial(poly, param, w)
        end
    end

    @testset "Example 4.1 (P3)" begin
        @var x y z
        sagbi = [x^2 + 1, y^2 + 1, x * y + z^2, 1]
        w = Int[-2, -2, -3]
        param = Parametrization(sagbi)

        @test degree_map(param; rng = MersenneTwister(1)) == 8
        @test degree_monomial_map(param, w; rng = MersenneTwister(1)) == 8

        @var p[1:4]
        got = nothing
        for seed in (20260210, 20260211, 20260212)
            rng = MersenneTwister(seed)
            lin_sys = [_nonzero_int_matrix_paper(rng, 3, 4) * p]
            linear = LinearSection(lin_sys)
            result = SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
            if length(result.solutions) == 8
                got = result
                break
            end
        end

        @test got !== nothing
        @test length(got.solutions) == 8
        @test length(solutions(got.tracker_result)) == 8
    end

    @testset "Example 4.2 (semimixed)" begin
        @var x y z
        sagbi = [[x, y, (x^2 + y^2), 1], [y, z, (x^2 + y^2), (x^3 + z^3)]]
        param = Parametrization(sagbi)
        w = detect_weight(param)

        @test degree_map(param; rng = MersenneTwister(5)) == 1
        @test degree_monomial_map(param, w; rng = MersenneTwister(5)) == 1

        rng = MersenneTwister(20260210)
        @var p[1:4]
        @var q[1:4]
        lin_sys = [
            _nonzero_int_matrix_paper(rng, 2, 4) * p,
            _nonzero_int_matrix_paper(rng, 1, 4) * q,
        ]

        linear = LinearSection(lin_sys)
        result = SagbiHomotopy.solve(linear, param; weight = w)
        @test length(result.solutions) == 6
        @test length(solutions(result.tracker_result)) == 6
    end

    @testset "Example 4.3 (base locus)" begin
        @var x y
        sagbi = [[x * (x^2 + y^2 - 2 * x), x * (5 - 4 * y), y * (x^2 + y^2 - 2 * x), y * (5 - 4 * y)]]
        w = Int[-2, -1]
        param = Parametrization(sagbi)

        @test degree_map(param; rng = MersenneTwister(7)) == 2
        @test degree_monomial_map(param, w; rng = MersenneTwister(7)) == 1

        rng = MersenneTwister(20260210)
        @var p[1:4]
        lin_sys = [rand(rng, ComplexF64, 2, 4) * p]
        linear = LinearSection(lin_sys)

        result = SagbiHomotopy.solve(linear, param; weight = w)
        @test length(result.solutions) == 2

        result2 = SagbiHomotopy.solve(linear, param; weight = w, include_base_locus = true)
        @test length(result2.solutions) == 4
        @test length(solutions(result2.tracker_result)) == 2
    end

    @testset "Example 3.5 (degree drop)" begin
        @var x y z
        sagbi = [x^3, y^3, x^6 * y^3 + x * z^2 + y^3]
        w = Int[-3, -3, -8]
        param = Parametrization(sagbi)

        @test degree_map(param; rng = MersenneTwister(11)) == 3
        @test degree_monomial_map(param, w; rng = MersenneTwister(11)) == 18
    end

    @testset "Section 5.1 (Gr(2,6) linear slice)" begin
        rng = MersenneTwister(20260210)
        gr = grassmannian_parametrization(2, 6)
        param = gr.parametrization
        w = gr.weight

        @var z[1:15]
        A = randn(rng, ComplexF64, 8, 15)
        linear = LinearSection(A * z)

        result = SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
        @test length(result.solutions) == 14
        @test length(solutions(result.tracker_result)) == 14
    end

    @testset "Section 5.1 (Gr(2,5) linear slice)" begin
        rng = MersenneTwister(20260210)
        gr = grassmannian_parametrization(2, 5)
        param = gr.parametrization
        w = gr.weight

        @var z[1:10]
        A = randn(rng, ComplexF64, 6, 10)
        linear = LinearSection(A * z)

        result = SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
        @test length(result.solutions) == 5
        @test length(solutions(result.tracker_result)) == 5
    end

    @testset "Section 5.1 (Gr(3,6) linear slice)" begin
        rng = MersenneTwister(20260210)
        gr = grassmannian_parametrization(3, 6)
        param = gr.parametrization
        w = gr.weight

        @var z[1:20]
        A = randn(rng, ComplexF64, 9, 20)
        linear = LinearSection(A * z)

        result = SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
        @test length(result.solutions) == 42
        @test length(solutions(result.tracker_result)) == 42
    end

    @testset "Section 5.2 (nonlinear resonator / oscillators, N=1 M=2)" begin
        rng = MersenneTwister(20260210)
        N = 1
        M = 2

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
            end
        end

        sagbi = [polynomials[i, :] for i in 1:(2 * M * N)]
        param = Parametrization(sagbi)

        c = rand(rng, Float64, 2 * M * N, L)
        lin_sys = [[sum(z[i, :] .* c[i, :])] for i in 1:(2 * M * N)]
        linear = LinearSection(lin_sys)

        w = vcat(fill(2, M), fill(1, M))
        result = SagbiHomotopy.solve(linear, param; weight = w, check_degree = false)
        @test length(result.solutions) == 25
        @test length(solutions(result.tracker_result)) == 25

        full_eqs = [subs(lin_sys[i][1], collect(z[i, :]) => sagbi[i]) for i in 1:length(sagbi)]
        res_full = HomotopyContinuation.solve(full_eqs; show_progress = false)
        @test length(HomotopyContinuation.solutions(res_full)) == 25
    end
end
