"""
    InputShapeError(message)

Exception thrown when user-facing container inputs violate the strict shape contract.
"""
struct InputShapeError <: Exception
    message::String
end

"""
    IncompatibleProblemError(message)

Exception thrown when linear-section and parametrization groups are structurally incompatible.
"""
struct IncompatibleProblemError <: Exception
    message::String
end

"""
    NoSagbiWeightError(message)

Exception thrown when no SAGBI-realizing weight vector can be detected.
"""
struct NoSagbiWeightError <: Exception
    message::String
end

"""
    DegreeDropError(degree_map, degree_monomial_map)

Exception thrown when degree diagnostics detect a drop and `allow_degree_drop=false`.
"""
struct DegreeDropError <: Exception
    degree_map::Int
    degree_monomial_map::Int
end

"""
    Base.showerror(io, err::InputShapeError)

Internal display helper for `InputShapeError`.
"""
function Base.showerror(io::IO, err::InputShapeError)
    return print(io, err.message)
end

"""
    Base.showerror(io, err::IncompatibleProblemError)

Internal display helper for `IncompatibleProblemError`.
"""
function Base.showerror(io::IO, err::IncompatibleProblemError)
    return print(io, err.message)
end

"""
    Base.showerror(io, err::NoSagbiWeightError)

Internal display helper for `NoSagbiWeightError`.
"""
function Base.showerror(io::IO, err::NoSagbiWeightError)
    return print(io, err.message)
end

"""
    Base.showerror(io, err::DegreeDropError)

Internal display helper for `DegreeDropError`.
"""
function Base.showerror(io::IO, err::DegreeDropError)
    return print(
        io,
        "degree of monomial parametrization drops from ",
        err.degree_map,
        " to ",
        err.degree_monomial_map,
        ". Set allow_degree_drop=true to continue.",
    )
end
