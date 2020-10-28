"""
    Geometry

This module defines the basic geometrical objects manipulated. 

Exports:
$(EXPORTS)
"""
module Geometry 

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype
import Base: union, setdiff, intersect, issubset

using DocStringExtensions
using StaticArrays
using LinearAlgebra

using WaveProp.Utils

include("point.jl")
include("hyperrectangle.jl")
include("domain.jl")
include("element.jl")
include("transformations.jl")
include("extrudedElement.jl")

export ElementaryEntity, Domain,
    dim, entities, tag, tags, remove, assertequaldim,
    boundary, skeleton, internal_boundary, external_boundary,
    GenericMesh, LagrangeElement, Point, ReferenceLine, ReferenceTriangle, ReferenceTetrahedron, ReferenceSquare,
    domain, jacobian, normal, LagrangeLine, LagrangeTriangle, LagrangeTetrahedron, AbstractElement, ambient_dimension, geometric_dimension, push_forward,
    IMT, Duffy, GeometricTransformation, bounding_box, HyperRectangle, center, diameter, radius, line, triangle, rectangle, extrude, translate, ParametricLine, AbstractReferenceShape, measure

end # module