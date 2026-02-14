"""
    _materialize_expression_vector(values, context)

Internal helper that validates `values` entries are symbolic (`Expression`/`Variable`) and
returns a concretely typed `Vector{Expression}`.
"""
function _materialize_expression_vector(values::AbstractVector, context::AbstractString)
    out = Vector{Expression}(undef, length(values))
    for i in eachindex(values)
        value = values[i]
        if !(value isa Expression || value isa Variable)
            throw(InputShapeError("$context must contain only HomotopyContinuation.Expression entries."))
        end
        out[i] = Expression(value)
    end
    return out
end

"""
    _normalize_groups(values, context)

Internal helper that enforces the strict one-level grouping contract used by public constructors.

Accepted shapes:
- `Vector{Expression}` (single group)
- `Vector{Vector{Expression}}` (multiple groups)

Anything deeper is rejected.
"""
function _normalize_groups(values::AbstractVector, context::AbstractString)
    isempty(values) && throw(InputShapeError("$context cannot be empty."))

    if all(v -> v isa Expression || v isa Variable, values)
        group = _materialize_expression_vector(values, context)
        isempty(group) && throw(InputShapeError("$context cannot be empty."))
        return [group]
    end

    groups = Vector{Vector{Expression}}(undef, length(values))
    for i in eachindex(values)
        group = values[i]
        if !(group isa AbstractVector)
            throw(InputShapeError("$context must be Vector{Expression} or Vector{Vector{Expression}}."))
        end
        if any(v -> v isa AbstractVector, group)
            throw(InputShapeError("$context does not allow nesting deeper than one group level."))
        end
        materialized = _materialize_expression_vector(group, context)
        isempty(materialized) && throw(InputShapeError("$context cannot contain empty groups."))
        groups[i] = materialized
    end

    return groups
end

"""
    _flatten_groups(groups)

Internal helper that concatenates grouped expressions into one flat vector.
"""
function _flatten_groups(groups::Vector{Vector{Expression}})
    total = sum(length, groups)
    out = Vector{Expression}(undef, total)
    idx = 1
    for group in groups
        for poly in group
            out[idx] = poly
            idx += 1
        end
    end
    return out
end

"""
    _validate_weight(weight, nvars)

Internal helper that checks weight length and returns a concrete `Vector{Int}`.
"""
function _validate_weight(weight::AbstractVector{<:Integer}, nvars::Int)
    if length(weight) != nvars
        throw(InputShapeError("weight length $(length(weight)) does not match number of variables $nvars."))
    end
    return Int[wi for wi in weight]
end
