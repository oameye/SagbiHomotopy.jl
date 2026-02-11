function _leading_index(weights::AbstractVector{<:Integer})
    lead = minimum(weights)
    idx = findfirst(==(lead), weights)
    return idx::Int, lead
end

function leading_monomial(poly::Expression, vars::Vector{Variable}, weight::AbstractVector{<:Integer})
    checked_weight = _validate_weight(weight, length(vars))
    exps, coeffs = exponents_coefficients(poly, vars)
    weighted = vec(transpose(exps) * checked_weight)
    idx, _ = _leading_index(weighted)

    term = coeffs[idx]
    for j in eachindex(vars)
        exp_ij = exps[j, idx]
        if exp_ij != 0
            term *= vars[j]^exp_ij
        end
    end

    return Expression(term)
end

function leading_monomial(poly::Expression, parametrization::Parametrization, weight::AbstractVector{<:Integer})
    return leading_monomial(poly, parametrization.vars, weight)
end

function _weight_deformation_with_t(
        poly::Expression,
        vars::Vector{Variable},
        weight::AbstractVector{<:Integer},
        t::Variable,
    )
    checked_weight = _validate_weight(weight, length(vars))

    exps, coeffs = exponents_coefficients(poly, vars)
    weighted = vec(transpose(exps) * checked_weight)
    _, lead = _leading_index(weighted)

    nterms = length(coeffs)
    first_term = coeffs[1]
    first_t_exp = weighted[1] - lead
    if first_t_exp != 0
        first_term *= t^first_t_exp
    end
    for j in eachindex(vars)
        exp_1j = exps[j, 1]
        if exp_1j != 0
            first_term *= vars[j]^exp_1j
        end
    end

    acc = first_term
    for i in 2:nterms
        term = coeffs[i]
        t_exp = weighted[i] - lead
        if t_exp != 0
            term *= t^t_exp
        end
        for j in eachindex(vars)
            exp_ij = exps[j, i]
            if exp_ij != 0
                term *= vars[j]^exp_ij
            end
        end
        acc += term
    end

    return Expression(acc)
end

function weight_deformation(
        poly::Expression,
        vars::Vector{Variable},
        weight::AbstractVector{<:Integer};
        t::Union{Nothing, Variable} = nothing,
    )
    if isnothing(t)
        @unique_var generated_t
        return _weight_deformation_with_t(poly, vars, weight, generated_t)
    end
    return _weight_deformation_with_t(poly, vars, weight, t)
end

function weight_deformation(
        poly::Expression,
        parametrization::Parametrization,
        weight::AbstractVector{<:Integer};
        t::Union{Nothing, Variable} = nothing,
    )
    return weight_deformation(poly, parametrization.vars, weight; t = t)
end

function degree_map(parametrization::Parametrization; rng::AbstractRNG = default_rng(), random_range::Int = 100)
    random_range > 0 || throw(InputShapeError("random_range must be positive."))

    ngroups = length(parametrization.groups)
    @unique_var s[1:ngroups]

    homogenized = Vector{Vector{Expression}}(undef, ngroups)
    for i in eachindex(parametrization.groups)
        group = parametrization.groups[i]
        homogenized_group = Vector{Expression}(undef, length(group))
        for j in eachindex(group)
            homogenized_group[j] = s[i] * group[j]
        end
        homogenized[i] = homogenized_group
    end

    hc_vars = Vector{Variable}(variables(_flatten_groups(homogenized)))
    ring, ring_vars_tuple = polynomial_ring(QQ, string.(hc_vars))
    ring_vars = Vector{QQMPolyRingElem}(ring_vars_tuple)

    mapped = Vector{QQMPolyRingElem}()
    for group in homogenized
        for poly in group
            push!(mapped, _hc_to_oscar(poly, ring_vars, hc_vars))
        end
    end

    values = Vector{typeof(QQ(1))}(undef, length(ring_vars))
    s_set = Set{Variable}(s)
    nonzero_nums = collect(-random_range:-1)
    append!(nonzero_nums, 1:random_range)

    for i in eachindex(hc_vars)
        if hc_vars[i] in s_set
            values[i] = QQ(1)
        else
            num = rand(rng, nonzero_nums)
            den = rand(rng, 1:random_range)
            values[i] = QQ(num // den)
        end
    end

    evaluated = Vector{QQMPolyRingElem}(undef, length(mapped))
    for i in eachindex(mapped)
        evaluated[i] = ring(Oscar.evaluate(mapped[i], ring_vars_tuple, values))
    end

    fiber_ideal = ideal(mapped .- evaluated)
    decomposition = absolute_primary_decomposition(fiber_ideal)

    return sum(component[4] for component in decomposition)
end

function degree_map(parametrization::AbstractVector; kwargs...)
    return degree_map(Parametrization(parametrization); kwargs...)
end

function degree_monomial_map(
        parametrization::Parametrization,
        weight::AbstractVector{<:Integer};
        rng::AbstractRNG = default_rng(),
        random_range::Int = 100,
    )
    checked_weight = _validate_weight(weight, length(parametrization.vars))

    monomial_groups = Vector{Vector{Expression}}(undef, length(parametrization.groups))
    for i in eachindex(parametrization.groups)
        group = parametrization.groups[i]
        monomial_group = Vector{Expression}(undef, length(group))
        for j in eachindex(group)
            monomial_group[j] = leading_monomial(group[j], parametrization.vars, checked_weight)
        end
        monomial_groups[i] = monomial_group
    end

    return degree_map(Parametrization(monomial_groups); rng = rng, random_range = random_range)
end

function degree_monomial_map(parametrization::AbstractVector, weight::AbstractVector{<:Integer}; kwargs...)
    return degree_monomial_map(Parametrization(parametrization), weight; kwargs...)
end

function _isolated_nsols_baselocus(sagbi_group::Vector{Expression})
    hc_vars = Vector{Variable}(variables(sagbi_group))
    _, ring_vars_tuple = polynomial_ring(QQ, string.(hc_vars))
    ring_vars = Vector{QQMPolyRingElem}(ring_vars_tuple)

    mapped = Vector{QQMPolyRingElem}(undef, length(sagbi_group))
    for i in eachindex(sagbi_group)
        mapped[i] = _hc_to_oscar(sagbi_group[i], ring_vars, hc_vars)
    end

    decomposition = absolute_primary_decomposition(ideal(mapped))

    isolated_degree = 0
    for component in decomposition
        if Oscar.dim(component[3]) == 0
            isolated_degree += component[4]
        end
    end

    return isolated_degree
end

function _base_locus_solutions(sagbi_group::Vector{Expression})
    return solutions(HomotopyContinuation.solve(sagbi_group; show_progress = false))
end
