function _initial_form(poly::QQMPolyRingElem, ring::QQMPolyRing, weight::AbstractVector{<:Integer})
    coeffs, exps = _coeffs_exponents(poly)
    isempty(coeffs) && return poly - poly

    ctx = MPolyBuildCtx(ring)
    push_term!(ctx, coeffs[1], exps[1])
    max_degree = _weighted_degree(exps[1], weight)

    for i in 2:length(coeffs)
        degree_i = _weighted_degree(exps[i], weight)
        if degree_i > max_degree
            finish(ctx)
            ctx = MPolyBuildCtx(ring)
            push_term!(ctx, coeffs[i], exps[i])
            max_degree = degree_i
        elseif degree_i == max_degree
            push_term!(ctx, coeffs[i], exps[i])
        end
    end

    return finish(ctx)
end

function _sagbi_criterion(G::Vector{QQMPolyRingElem}, ring::QQMPolyRing, weight::Vector{Int})
    n = length(G)
    n == 0 && return false

    order = wdeglex(ring, weight)
    leading_terms = Vector{QQMPolyRingElem}(undef, n)
    for i in eachindex(G)
        leading_terms[i] = Oscar.leading_term(G[i], ordering = order)
    end

    graded_ring, _ = Oscar.graded_polynomial_ring(Oscar.QQ, "z" .* string.(1:n))
    map_leading = hom(graded_ring, ring, leading_terms)
    map_full = hom(graded_ring, ring, G)

    ideal_leading = kernel(map_leading)
    ideal_full = kernel(map_full)

    quotient_leading, _ = quo(graded_ring, ideal_leading)
    quotient_full, _ = quo(graded_ring, ideal_full)

    return string(Oscar.hilbert_series(quotient_leading)) == string(Oscar.hilbert_series(quotient_full))
end

function _weight_vector_realizing_sagbi(G::Vector{QQMPolyRingElem}, ring::QQMPolyRing)
    isempty(G) && return nothing

    newton_poly = Oscar.newton_polytope(reduce(*, G; init = one(ring)))
    vertices = Oscar.vertices(newton_poly)
    isempty(vertices) && return nothing

    normal_cones = maximal_cones(normal_fan(newton_poly))

    n = length(vertices[1])
    identity = Oscar.identity_matrix(Oscar.ZZ, n)
    nonnegative_orthant = Oscar.positive_hull(-identity)

    for cone in normal_cones
        intersection = Oscar.intersect(cone, nonnegative_orthant)
        if Oscar.dim(intersection) > 0
            rays = matrix(Oscar.ZZ, Oscar.rays(intersection))
            weight = vec(ones(Int, size(rays, 1)) * rays)

            if all(wi -> wi < 0, weight)
                candidate = Int[wi for wi in weight]
                if _sagbi_criterion(G, ring, -candidate)
                    return candidate
                end
            end
        end
    end

    return nothing
end

function detect_weight(parametrization::Parametrization)
    ngroups = length(parametrization.groups)
    vars = parametrization.vars

    @unique_var t[1:ngroups]
    names = vcat(["t$i" for i in 1:ngroups], string.(vars))
    ring, ring_vars = Oscar.polynomial_ring(Oscar.QQ, names)

    t_vars = ring_vars[1:ngroups]
    x_vars = ring_vars[(ngroups + 1):end]

    mapped = Vector{QQMPolyRingElem}()
    for i in eachindex(parametrization.groups)
        for poly in parametrization.groups[i]
            push!(mapped, t_vars[i] * _hc_to_oscar(poly, x_vars, vars))
        end
    end

    weight = _weight_vector_realizing_sagbi(mapped, ring)
    if isnothing(weight)
        throw(NoSagbiWeightError("No weight realizes the provided parametrization as a SAGBI basis."))
    end

    return weight[(ngroups + 1):end]
end

function detect_weight(parametrization::AbstractVector)
    return detect_weight(Parametrization(parametrization))
end
