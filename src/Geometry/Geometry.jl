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

using WaveProp.Utils

include("point.jl")
include("domain.jl")
include("element.jl")
include("transformations.jl")
include("mesh.jl")

export ElementaryEntity, Domain,
    dim, entities, tag, tags, remove, assertequaldim,
    boundary, skeleton, internal_boundary, external_boundary,
    GenericMesh, LagrangeElement, Point, ReferenceLine, ReferenceTriangle, ReferenceTetrahedron, ReferenceSquare,
    reference_element, jacobian, LagrangeLine, LagrangeTriangle, LagrangeTetrahedron, etypes, AbstractElement, ambient_dimension, geometric_dimension, 
    IMT, Duffy, GeometricTransformation

end # module