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

using WaveProp
using WaveProp.Utils

WaveProp.@import_interface

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
    Domain,
    ReferencePoint,
    ReferenceLine,
    ReferenceTriangle,
    ReferenceTetrahedron,
    ReferenceSquare,
    HyperRectangle,
    # functions
    clear_entities!,
    tag,
    key,
    assertequaldim,
    boundary,
    skeleton,
    internal_boundary,
    external_boundary,
    flip_normal,
    low_corner,
    high_corner,
    measure,
    vertices,
    line,
    # global variables
    TAGS,
    ENTITIES

end # module
