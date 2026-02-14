module SagbiHomotopy

using LinearAlgebra: I, det, dot, norm
using AbstractAlgebra: hom, ideal, kernel, matrix, polynomial_ring, quo
using AbstractAlgebra.Generic: MPolyBuildCtx, finish, push_term!

using Combinatorics: combinations

using HomotopyContinuation: HomotopyContinuation, evaluate, solutions
using HomotopyContinuation.ModelKit: @unique_var, @var, Expression, System, Variable, exponents_coefficients, support_coefficients

using MultivariatePolynomials: subs, variables
using Nemo: QQ, QQMPolyRing, QQMPolyRingElem

using Oscar: Oscar, absolute_primary_decomposition, maximal_cones, normal_fan
using Oscar.Orderings: wdeglex

using Random: Random, AbstractRNG, default_rng

export Parametrization,
    LinearSection,
    SagbiProblem,
    SolveOptions,
    SolveResult,
    SolveMetadata,
    InputShapeError,
    IncompatibleProblemError,
    NoSagbiWeightError,
    DegreeDropError,
    solve,
    detect_weight,
    possible_sagbi,
    degree_map,
    degree_monomial_map,
    leading_monomial,
    weight_deformation,
    grassmannian_parametrization

include("errors.jl")
include("normalize.jl")
include("types.jl")
include("convert.jl")
include("detection.jl")
include("degree.jl")
include("solve.jl")
include("convenience.jl")

end # module SagbiHomotopy
