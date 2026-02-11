struct InputShapeError <: Exception
    message::String
end

struct IncompatibleProblemError <: Exception
    message::String
end

struct NoSagbiWeightError <: Exception
    message::String
end

struct DegreeDropError <: Exception
    degree_map::Int
    degree_monomial_map::Int
end

function Base.showerror(io::IO, err::InputShapeError)
    return print(io, err.message)
end

function Base.showerror(io::IO, err::IncompatibleProblemError)
    return print(io, err.message)
end

function Base.showerror(io::IO, err::NoSagbiWeightError)
    return print(io, err.message)
end

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
