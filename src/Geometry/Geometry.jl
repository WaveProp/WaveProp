"""
    module Geometry

Module defines basic geometrical concepts.
"""
module Geometry

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype, keys
import Base: union, setdiff, intersect, issubset

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements
using OrderedCollections

# import all methods in WaveProp.INTERFACE_METHODS
using WaveProp
WaveProp.@import_interface

using WaveProp.Utils

include("hyperrectangle.jl")
include("referenceshapes.jl")
include("entities.jl")
include("domain.jl")
include("simpleshapes.jl")

export
    # abstract types
    AbstractReferenceShape,
    AbstractEntity,
    AbstractParametricBody,
    # types
    ElementaryEntity,
    ParametricEntity,
    ParametricCurve,
    ClosedEntity,
    Circle,
    Kite,
    Domain,
    ReferenceLine,
    ReferenceTriangle,
    ReferenceTetrahedron,
    ReferenceSquare,
    HyperRectangle,
    SVector,
    # functions
    clear_entities!,
    tag,
    key,
    remove,
    assertequaldim,
    boundary,
    skeleton,
    internal_boundary,
    external_boundary,
    domain,
    parametrization,
    jacobian,
    normal,
    flip_normal,
    ambient_dimension,
    geometric_dimension,
    bounding_box,
    center,
    diameter,
    radius,
    low_corner,
    high_corner,
    measure,
    number_of_nodes,
    vertices,
    line,
    # global variables
    TAGS,
    ENTITIES

end # module
