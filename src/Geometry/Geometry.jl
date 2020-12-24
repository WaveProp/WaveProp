"""
    module Geometry

This module defines basic geometrical objects.
"""
module Geometry 

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype, keys
import Base: union, setdiff, intersect, issubset

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements

using WaveProp.Utils

include("point.jl")
include("hyperrectangle.jl")
include("referenceshapes.jl")
include("entities.jl")
include("domain.jl")
include("parametricbody.jl")

export 
    # abstract types
    AbstractReferenceShape,
    AbstractEntity,
    AbstractParametricBody, 
    # types
    ElementaryEntity, 
    ParametricEntity,
    ParametricBody,
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
    entities, 
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
    ambient_dimension, 
    geometric_dimension, 
    bounding_box, 
    center, 
    diameter, 
    radius, 
    low_corner,
    high_corner,
    # global variables
    TAGS,
    ENTITIES
        
end # module
