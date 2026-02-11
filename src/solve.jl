function _solution_signature(sol, digits::Int)
    coords = collect(sol)
    io = IOBuffer()
    for i in eachindex(coords)
        z = coords[i]
        if i > 1
            print(io, ';')
        end
        print(io, round(real(z); digits = digits), ',', round(imag(z); digits = digits))
    end
    return String(take!(io))
end

function _build_deformed_groups(parametrization::Parametrization, weight::Vector{Int}, t::Variable)
    groups = Vector{Vector{Expression}}(undef, length(parametrization.groups))
    for i in eachindex(parametrization.groups)
        group = parametrization.groups[i]
        deformed = Vector{Expression}(undef, length(group))
        for j in eachindex(group)
            deformed[j] = weight_deformation(group[j], parametrization.vars, weight; t = t)
        end
        groups[i] = deformed
    end
    return groups
end

function _append_substituted_equations!(
        equations::Vector{Expression},
        linear_group::Vector{Expression},
        linear_vars::Vector{Variable},
        deformed_group::Vector{Expression},
    )
    substituted = subs(linear_group, linear_vars => deformed_group)
    append!(equations, substituted)
    return nothing
end

function _augment_base_locus!(
        solutions_vec::Vector,
        problem::SagbiProblem,
        system::System,
        opts::SolveOptions,
    )
    signatures = Set{String}()
    for sol in solutions_vec
        push!(signatures, _solution_signature(sol, opts.signature_digits))
    end

    candidates = 0
    kept = 0

    for group in problem.parametrization.groups
        n_isolated = _isolated_nsols_baselocus(group)
        if n_isolated == 0
            continue
        end

        missing_solutions = _base_locus_solutions(group)
        for missing_sol in missing_solutions
            candidates += 1
            coords = collect(missing_sol)

            if any(z -> abs(z) < opts.atol, coords)
                continue
            end

            residual = norm(evaluate(system, missing_sol, [1]), Inf)
            if residual >= opts.atol
                continue
            end

            sig = _solution_signature(missing_sol, opts.signature_digits)
            if sig in signatures
                continue
            end

            push!(signatures, sig)
            push!(solutions_vec, missing_sol)
            kept += 1
        end
    end

    return candidates, kept
end

function solve(problem::SagbiProblem, opts::SolveOptions = SolveOptions())
    opts.random_range > 0 || throw(InputShapeError("random_range must be positive."))

    weight = if isnothing(opts.weight)
        detect_weight(problem.parametrization)
    else
        _validate_weight(opts.weight, length(problem.parametrization.vars))
    end

    degree_full = nothing
    degree_monomial = nothing
    degree_drop = false

    if opts.check_degree
        degree_full = degree_map(problem.parametrization; rng = opts.rng, random_range = opts.random_range)
        degree_monomial = degree_monomial_map(
            problem.parametrization,
            weight;
            rng = opts.rng,
            random_range = opts.random_range,
        )
        degree_drop = degree_full > degree_monomial

        if degree_drop && !opts.allow_degree_drop
            throw(DegreeDropError(degree_full, degree_monomial))
        end
    end

    @unique_var t
    deformed_groups = _build_deformed_groups(problem.parametrization, weight, t)

    total_equations = sum(length, problem.linear.groups)
    equations = Vector{Expression}()
    sizehint!(equations, total_equations)

    scale = randn(opts.rng, ComplexF64)

    for i in eachindex(problem.linear.groups)
        linear_group = problem.linear.groups[i]
        linear_vars = problem.linear.vars_per_group[i]

        if opts.vary_linear_part
            nrows = length(linear_group)
            nvars = length(linear_vars)
            random_matrix = randn(opts.rng, ComplexF64, nrows, nvars)
            random_linear = random_matrix * linear_vars

            varied_group = Vector{Expression}(undef, nrows)
            for j in eachindex(varied_group)
                varied_group[j] = (1 - t) * random_linear[j] + (scale * t) * linear_group[j]
            end

            _append_substituted_equations!(equations, varied_group, linear_vars, deformed_groups[i])
        else
            _append_substituted_equations!(equations, linear_group, linear_vars, deformed_groups[i])
        end
    end

    system = System(equations, variables = problem.parametrization.vars, parameters = [t])

    start_result = HomotopyContinuation.solve(system; target_parameters = [0], show_progress = false)
    start_solutions = solutions(start_result)

    tracker_result = HomotopyContinuation.solve(
        system,
        start_solutions;
        start_parameters = [0],
        target_parameters = [1],
        show_progress = false,
    )

    solutions_vec = collect(solutions(tracker_result))

    base_candidates = 0
    base_kept = 0
    if opts.include_base_locus
        base_candidates, base_kept = _augment_base_locus!(solutions_vec, problem, system, opts)
    end

    metadata = SolveMetadata(weight, degree_full, degree_monomial, degree_drop, base_candidates, base_kept)

    return SolveResult(tracker_result, solutions_vec, metadata)
end

function solve(linear::LinearSection, parametrization::Parametrization; kwargs...)
    return solve(SagbiProblem(linear, parametrization), SolveOptions(; kwargs...))
end

function solve(linear::LinearSection, parametrization::Parametrization, opts::SolveOptions)
    return solve(SagbiProblem(linear, parametrization), opts)
end

function solve(linear::AbstractVector, parametrization::AbstractVector; kwargs...)
    return solve(SagbiProblem(linear, parametrization), SolveOptions(; kwargs...))
end
