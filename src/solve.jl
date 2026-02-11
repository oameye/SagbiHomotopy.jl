"""
    _solution_signature(sol, digits::Int)

Internal helper for robust numerical deduplication of solutions.

Builds a deterministic string signature by rounding real/imaginary parts of each coordinate
to `digits` decimal places.
"""
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

"""
    _build_deformed_groups(parametrization::Parametrization, weight::Vector{Int}, t::Variable)

Internal helper that applies [`weight_deformation`](@ref) to every polynomial in every
parametrization group, using homotopy parameter `t`.
"""
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

"""
    _append_substituted_equations!(equations, linear_group, linear_vars, deformed_group)

Internal helper that substitutes `linear_vars => deformed_group` into `linear_group` and appends
all resulting equations into `equations`.
"""
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

"""
    _augment_base_locus!(solutions_vec, problem, system, opts)

Internal helper for optional base-locus augmentation.

It computes base-locus solutions per parametrization group, filters them with torus and residual
checks, deduplicates numerically, appends accepted solutions to `solutions_vec`, and returns
`(candidates, kept)` counters.
"""
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

        missing_solutions = _base_locus_solutions(group; show_progress = opts.show_progress)
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

"""
    solve(problem::SagbiProblem, opts::SolveOptions=SolveOptions())

Solve a validated SAGBI homotopy problem and return a [`SolveResult`](@ref).

This is the main entrypoint. It composes grouped linear equations with a grouped
parametrization, constructs a one-parameter homotopy based on a SAGBI-induced weight
deformation, tracks solutions from `t=0` to `t=1`, and optionally augments with base-locus
solutions.

The method enforces strict typing and deterministic control through [`SolveOptions`](@ref).
If no weight is provided, [`detect_weight`](@ref) is called.

# Workflow
1. Validate runtime options (`random_range > 0`)
2. Determine effective weight (`opts.weight` or `detect_weight`)
3. Optionally run degree diagnostics (`check_degree`)
4. Build deformed parametrization and substituted system
5. Solve start system at `t=0`
6. Track to `t=1`
7. Optionally augment with base-locus points (`include_base_locus`)
8. Return `SolveResult(tracker_result, solutions, metadata)`

# Arguments
- `problem::SagbiProblem`: validated linear-section + parametrization bundle.
- `opts::SolveOptions`: solver and diagnostic options.

# Keyword behavior (via `SolveOptions`)
- `weight`: fixed integer weight vector; if omitted, autodetected.
- `check_degree`: compute map degree and monomial map degree.
- `allow_degree_drop`: permit continuation when degree drops.
- `include_base_locus`: append validated base-locus solutions.
- `vary_linear_part`: use randomized linear interpolation in the homotopy.
- `rng`: random number generator for all randomized branches.
- `random_range`: integer sampling range for degree estimation.
- `atol`: numerical tolerance used in base-locus acceptance checks.
- `signature_digits`: dedup precision for numerical solution signatures.
- `show_progress`: controls progress output in HomotopyContinuation solves.

# Returns
- [`SolveResult`](@ref): contains raw tracker result, final solution vector, and metadata.

# Errors
- [`InputShapeError`](@ref): invalid runtime option values.
- [`NoSagbiWeightError`](@ref): no SAGBI-realizing weight found (when weight autodetection fails).
- [`DegreeDropError`](@ref): degree drop detected with `check_degree=true` and
  `allow_degree_drop=false`.
"""
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

    start_result = HomotopyContinuation.solve(system; target_parameters = [0], show_progress = opts.show_progress)
    start_solutions = solutions(start_result)

    tracker_result = HomotopyContinuation.solve(
        system,
        start_solutions;
        start_parameters = [0],
        target_parameters = [1],
        show_progress = opts.show_progress,
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

"""
    solve(linear::LinearSection, parametrization::Parametrization; kwargs...)

Convenience wrapper around [`solve(problem::SagbiProblem, opts::SolveOptions)`](@ref).

Builds a [`SagbiProblem`](@ref) from typed inputs and a [`SolveOptions`](@ref) from `kwargs`.
"""
function solve(linear::LinearSection, parametrization::Parametrization; kwargs...)
    return solve(SagbiProblem(linear, parametrization), SolveOptions(; kwargs...))
end

"""
    solve(linear::LinearSection, parametrization::Parametrization, opts::SolveOptions)

Convenience wrapper that solves typed inputs with an explicit `SolveOptions` instance.
"""
function solve(linear::LinearSection, parametrization::Parametrization, opts::SolveOptions)
    return solve(SagbiProblem(linear, parametrization), opts)
end

"""
    solve(linear::AbstractVector, parametrization::AbstractVector; kwargs...)

High-level convenience wrapper.

This method first constructs [`LinearSection`](@ref) and [`Parametrization`](@ref), validates
compatibility via [`SagbiProblem`](@ref), then dispatches to [`solve`](@ref).
"""
function solve(linear::AbstractVector, parametrization::AbstractVector; kwargs...)
    return solve(SagbiProblem(linear, parametrization), SolveOptions(; kwargs...))
end
