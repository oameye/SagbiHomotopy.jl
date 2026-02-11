"""
    Parametrization{E<:Expression}

Typed container for a grouped SAGBI parametrization.

# Fields
- `groups::Vector{Vector{E}}`: parametrization groups. Each inner vector corresponds to one
  parameter block and contains symbolic expressions used as substitution targets.
- `vars::Vector{Variable}`: cached global variable list appearing in `groups`, used to avoid
  repeated variable discovery in downstream algorithms.
"""
struct Parametrization{E <: Expression}
    groups::Vector{Vector{E}}
    vars::Vector{Variable}
end

"""
    LinearSection{E<:Expression}

Typed container for grouped linear equations to be composed with a parametrization.

# Fields
- `groups::Vector{Vector{E}}`: grouped linear equations.
- `vars_per_group::Vector{Vector{Variable}}`: cached variables per group, used to validate
  compatibility with matching parametrization groups.
"""
struct LinearSection{E <: Expression}
    groups::Vector{Vector{E}}
    vars_per_group::Vector{Vector{Variable}}
end

"""
    SagbiProblem{E<:Expression}

Validated problem bundle combining a [`LinearSection`](@ref) and [`Parametrization`](@ref).

A `SagbiProblem` guarantees:
- same number of groups in `linear` and `parametrization`
- for each group `i`, `length(linear.vars_per_group[i]) == length(parametrization.groups[i])`

Those checks are enforced by the inner constructor.
"""
struct SagbiProblem{E <: Expression}
    linear::LinearSection{E}
    parametrization::Parametrization{E}

    function SagbiProblem{E}(linear::LinearSection{E}, parametrization::Parametrization{E}) where {E <: Expression}
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

"""
    SolveOptions

Configuration for [`solve`](@ref).

# Fields
- `weight::Union{Nothing,Vector{Int}}=nothing`: explicit SAGBI weight. If `nothing`,
  [`detect_weight`](@ref) is run.
- `check_degree::Bool=false`: whether to compute and compare `degree_map` and
  `degree_monomial_map`.
- `allow_degree_drop::Bool=false`: when `check_degree=true`, controls whether a detected
  degree drop throws [`DegreeDropError`](@ref).
- `include_base_locus::Bool=false`: whether to augment tracker solutions with validated
  base-locus candidates.
- `vary_linear_part::Bool=false`: whether to randomize linear equations along the homotopy.
- `rng::AbstractRNG=Random.default_rng()`: RNG used by all randomized branches.
- `random_range::Int=100`: sampling range used by degree estimation.
- `atol::Float64=1e-10`: tolerance for base-locus filtering and residual checks.
- `signature_digits::Int=8`: rounding precision used for numerical solution deduplication.
"""
Base.@kwdef struct SolveOptions{R <: AbstractRNG}
    weight::Union{Nothing, Vector{Int}} = nothing
    check_degree::Bool = false
    allow_degree_drop::Bool = false
    include_base_locus::Bool = false
    vary_linear_part::Bool = false
    rng::R = default_rng()
    random_range::Int = 100
    atol::Float64 = 1.0e-10
    signature_digits::Int = 8
end

"""
    SolveMetadata

Summary diagnostics returned by [`solve`](@ref).

# Fields
- `weight`: effective weight used for the run.
- `degree_map`: degree of the original map (or `nothing` when `check_degree=false`).
- `degree_monomial_map`: degree of the monomialized map (or `nothing`).
- `degree_drop_detected`: whether degree drop was detected.
- `base_locus_candidates`: number of base-locus candidates examined.
- `base_locus_kept`: number of candidates accepted and appended.
"""
struct SolveMetadata
    weight::Vector{Int}
    degree_map::Union{Nothing, Int}
    degree_monomial_map::Union{Nothing, Int}
    degree_drop_detected::Bool
    base_locus_candidates::Int
    base_locus_kept::Int
end

"""
    SolveResult{R,S}

Result object returned by [`solve`](@ref).

# Fields
- `tracker_result`: raw `HomotopyContinuation.solve` result.
- `solutions::Vector{S}`: final solution list (possibly base-locus augmented).
- `metadata::SolveMetadata`: diagnostic metadata for the run.
"""
struct SolveResult{R, S}
    tracker_result::R
    solutions::Vector{S}
    metadata::SolveMetadata
end

"""
    Parametrization(groups::AbstractVector)

Construct a validated [`Parametrization`](@ref) from grouped symbolic expressions.

Accepted input shapes are strict:
- `Vector{Expression}` (single group)
- `Vector{<:Vector{Expression}}` (multiple groups)

`Variable` entries are accepted and converted to `Expression`.

# Errors
- Throws [`InputShapeError`](@ref) when the shape is empty, nested deeper than one level,
  or contains non-symbolic entries.
"""
function Parametrization(groups::AbstractVector)
    normalized = _normalize_groups(groups, "Parametrization")
    vars_cache = Vector{Variable}(variables(_flatten_groups(normalized)))
    return Parametrization{Expression}(normalized, vars_cache)
end

"""
    LinearSection(groups::AbstractVector)

Construct a validated [`LinearSection`](@ref) from grouped linear symbolic equations.

Accepted input shapes are strict:
- `Vector{Expression}` (single group)
- `Vector{<:Vector{Expression}}` (multiple groups)

`Variable` entries are accepted and converted to `Expression`.

# Errors
- Throws [`InputShapeError`](@ref) for invalid shape/content (same rules as
  [`Parametrization`](@ref)).
"""
function LinearSection(groups::AbstractVector)
    normalized = _normalize_groups(groups, "LinearSection")
    vars_per_group = Vector{Vector{Variable}}(undef, length(normalized))
    for i in eachindex(normalized)
        vars_per_group[i] = Vector{Variable}(variables(normalized[i]))
    end
    return LinearSection{Expression}(normalized, vars_per_group)
end

"""
    SagbiProblem(linear::LinearSection, parametrization::Parametrization)

Create a validated [`SagbiProblem`](@ref) from typed inputs.

Validation enforces group-count and per-group arity compatibility.

# Errors
- Throws [`IncompatibleProblemError`](@ref) if groups are incompatible.
"""
function SagbiProblem(linear::LinearSection{E}, parametrization::Parametrization{E}) where {E <: Expression}
    return SagbiProblem{E}(linear, parametrization)
end

"""
    SagbiProblem(linear::AbstractVector, parametrization::AbstractVector)

Convenience constructor that first builds [`LinearSection`](@ref) and [`Parametrization`](@ref),
then validates compatibility.
"""
function SagbiProblem(linear::AbstractVector, parametrization::AbstractVector)
    return SagbiProblem(LinearSection(linear), Parametrization(parametrization))
end
