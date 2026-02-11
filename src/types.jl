struct Parametrization{E<:Expression}
    groups::Vector{Vector{E}}
    vars::Vector{Variable}
end

struct LinearSection{E<:Expression}
    groups::Vector{Vector{E}}
    vars_per_group::Vector{Vector{Variable}}
end

struct SagbiProblem{E<:Expression}
    linear::LinearSection{E}
    parametrization::Parametrization{E}

    function SagbiProblem{E}(linear::LinearSection{E}, parametrization::Parametrization{E}) where {E<:Expression}
        if length(linear.groups) != length(parametrization.groups)
            throw(IncompatibleProblemError("LinearSection and Parametrization must have the same number of groups."))
        end

        for i in eachindex(linear.groups)
            nvars_linear = length(linear.vars_per_group[i])
            nparams = length(parametrization.groups[i])
            if nvars_linear != nparams
                throw(
                    IncompatibleProblemError(
                        "Group $i mismatch: linear group has $nvars_linear variables, parametrization group has $nparams entries.",
                    ),
                )
            end
        end

        return new{E}(linear, parametrization)
    end
end

Base.@kwdef struct SolveOptions{R<:AbstractRNG}
    weight::Union{Nothing,Vector{Int}} = nothing
    check_degree::Bool = false
    allow_degree_drop::Bool = false
    include_base_locus::Bool = false
    vary_linear_part::Bool = false
    rng::R = default_rng()
    random_range::Int = 100
    atol::Float64 = 1.0e-10
    signature_digits::Int = 8
end

struct SolveMetadata
    weight::Vector{Int}
    degree_map::Union{Nothing,Int}
    degree_monomial_map::Union{Nothing,Int}
    degree_drop_detected::Bool
    base_locus_candidates::Int
    base_locus_kept::Int
end

struct SolveResult{R,S}
    tracker_result::R
    solutions::Vector{S}
    metadata::SolveMetadata
end

function Parametrization(groups::AbstractVector)
    normalized = _normalize_groups(groups, "Parametrization")
    vars_cache = Vector{Variable}(variables(_flatten_groups(normalized)))
    return Parametrization{Expression}(normalized, vars_cache)
end

function LinearSection(groups::AbstractVector)
    normalized = _normalize_groups(groups, "LinearSection")
    vars_per_group = Vector{Vector{Variable}}(undef, length(normalized))
    for i in eachindex(normalized)
        vars_per_group[i] = Vector{Variable}(variables(normalized[i]))
    end
    return LinearSection{Expression}(normalized, vars_per_group)
end

function SagbiProblem(linear::LinearSection{E}, parametrization::Parametrization{E}) where {E<:Expression}
    return SagbiProblem{E}(linear, parametrization)
end

function SagbiProblem(linear::AbstractVector, parametrization::AbstractVector)
    return SagbiProblem(LinearSection(linear), Parametrization(parametrization))
end
