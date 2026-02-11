function _monomial_from_support(vars_sys::Vector{Variable}, support_matrix::AbstractMatrix, col::Int)
    nvars = length(vars_sys)
    term = vars_sys[1]^0
    for j in 1:nvars
        exp_j = support_matrix[j, col]
        if exp_j != 0
            term *= vars_sys[j]^exp_j
        end
    end
    return term
end

function possible_sagbi(sys::System)
    supports, coeff_vectors = support_coefficients(sys)
    vars_sys = Vector{Variable}(variables(sys))

    linear_equations = Vector{Expression}(undef, length(supports))
    all_sagbi_terms = Vector{Expression}()
    seen_terms = Dict{Expression, Int}()

    for i in eachindex(supports)
        support_matrix = supports[i]
        coeffs = coeff_vectors[i]

        unique_coeffs = Vector{typeof(coeffs[1])}()
        for coeff in coeffs
            if !(coeff in unique_coeffs)
                push!(unique_coeffs, coeff)
            end
        end

        sagbi_part = Vector{Expression}(undef, length(unique_coeffs))
        for j in eachindex(unique_coeffs)
            coeff = unique_coeffs[j]
            idxs = findall(==(coeff), coeffs)

            first_idx = idxs[1]
            acc = _monomial_from_support(vars_sys, support_matrix, first_idx)
            for idx in idxs[2:end]
                acc += _monomial_from_support(vars_sys, support_matrix, idx)
            end

            sagbi_part[j] = acc
            if !haskey(seen_terms, acc)
                push!(all_sagbi_terms, acc)
                seen_terms[acc] = length(all_sagbi_terms)
            end
        end

        lin_eq = unique_coeffs[1] * sagbi_part[1]
        for j in 2:length(unique_coeffs)
            lin_eq += unique_coeffs[j] * sagbi_part[j]
        end

        linear_equations[i] = lin_eq
    end

    k = length(all_sagbi_terms)
    @var z[1:k]
    linear_in_new_vars = subs(linear_equations, all_sagbi_terms => collect(z))

    return LinearSection(linear_in_new_vars), Parametrization(all_sagbi_terms)
end

function grassmannian_parametrization(k::Int, m::Int)
    if !(1 <= k <= m)
        throw(InputShapeError("Expected 1 <= k <= m, got k=$k and m=$m."))
    end

    @var x[1:(k * (m - k))]
    matrix = hcat(Matrix(I, k, k), transpose(reshape(x, m - k, k)))

    index_sets = collect(combinations(1:m, k))
    sagbi = Vector{Expression}(undef, length(index_sets))
    for i in eachindex(index_sets)
        sagbi[i] = det(matrix[:, index_sets[i]])
    end

    weight = Int[i * (m - j + 1) for i in 0:(k - 1) for j in 1:(m - k)]

    return (
        vars = Vector{Variable}(x),
        parametrization = Parametrization(sagbi),
        weight = weight,
    )
end
