function _weighted_degree(monom_exponent::AbstractVector{<:Integer}, weight::AbstractVector{<:Integer})
    return dot(monom_exponent, weight)
end

function _coeffs_exponents(poly::QQMPolyRingElem)
    coeffs = collect(Oscar.coefficients(poly))
    exps = collect(Oscar.exponents(poly))
    return coeffs, exps
end

function _hc_to_oscar(
        poly::Expression,
        oscar_vars::AbstractVector{<:QQMPolyRingElem},
        hc_vars::AbstractVector{<:Variable},
    )
    exps, coeffs = exponents_coefficients(poly, hc_vars)

    zero_poly = oscar_vars[1] - oscar_vars[1]
    acc = zero_poly

    for i in eachindex(coeffs)
        term = coeffs[i]
        for j in eachindex(oscar_vars)
            exp_ij = exps[j, i]
            if exp_ij != 0
                term *= oscar_vars[j]^exp_ij
            end
        end
        acc += term
    end

    return acc
end
