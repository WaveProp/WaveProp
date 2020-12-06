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

"""
    const TAGS::Dict{Int,Vector{Int}}

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""    
const TAGS = Dict{Int,Vector{Int}}()

function _new_tag(dim)
    if !haskey(TAGS,dim) 
        return 1
    else
        tnew = maximum(TAGS[dim]) + 1
        return tnew    
    end
end

function _add_tag!(dim,tag)
    if is_new_tag(dim,tag) 
        if haskey(TAGS,dim)   
            push!(TAGS[dim],tag)
        else
            TAGS[dim] = [tag,]
        end
    end
    return TAGS    
end 

function is_new_tag(dim,tag)
    if haskey(TAGS,dim)
        existing_tags = TAGS[dim]
        if in(tag,existing_tags) 
            msg  = "entity of dimension $dim and tag $tag already exists in TAGS. 
            Creating a possibly duplicate entity."    
            @debug msg
            return false
        end
    end
    return true
end  

clear_tags!() = empty!(TAGS)

export ElementaryEntity, Domain,
    dim, entities, tag, tags, remove, assertequaldim,
    boundary, skeleton, internal_boundary, external_boundary,
    GenericMesh, LagrangeElement, Point, ReferenceLine, ReferenceTriangle, ReferenceTetrahedron, ReferenceSquare,
    domain, jacobian, normal, LagrangeLine, LagrangeTriangle, LagrangeTetrahedron, AbstractElement, ambient_dimension, geometric_dimension, push_forward,
    IMT, Duffy, GeometricTransformation, bounding_box, HyperRectangle, center, diameter, radius, line, triangle, rectangle, extrude, translate, ParametricLine, AbstractReferenceShape, measure,
    ParametricBody, AbstractParametricBody, ParametricEntity, SVector, SMatrix, Circle, Kite,AbstractParametricBody
end # module
