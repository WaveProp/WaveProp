"""
    Geometry

This module defines the basic geometrical objects manipulated in `WaveProp`.
"""
module Geometry 

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype
import Base: union, setdiff, intersect, issubset

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements

using WaveProp.Utils

include("point.jl")
include("hyperrectangle.jl")
include("domain.jl")
include("element.jl")
include("parametricentity.jl")
include("parametricbody.jl")

export ElementaryEntity, Domain,
    dim, entities, tag, tags, remove, assertequaldim,
    boundary, skeleton, internal_boundary, external_boundary,
    GenericMesh, LagrangeElement, Point, ReferenceLine, ReferenceTriangle, ReferenceTetrahedron, ReferenceSquare,
    domain, jacobian, normal, LagrangeLine, LagrangeTriangle, LagrangeTetrahedron, AbstractElement, ambient_dimension, geometric_dimension, push_forward,
    IMT, Duffy, GeometricTransformation, bounding_box, HyperRectangle, center, diameter, radius, line, triangle, rectangle, extrude, translate, ParametricLine, AbstractReferenceShape, measure,
    ParametricBody, AbstractParametricBody, ParametricEntity, SVector, SMatrix
end # module
