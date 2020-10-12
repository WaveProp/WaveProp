"""
    Geometry

This module defines the basic geometrical objects manipulated. 

Exports:
$(EXPORTS)
"""
module Geometry 

import Base: ==, length, in, iterate, getindex, lastindex, isempty
import Base: union, setdiff, intersect, issubset

using DocStringExtensions
using StaticArrays

using WaveProp.Utils

include("point.jl")
include("domain.jl")
include("element.jl")
include("mesh.jl")

export ElementaryEntity, Domain,
    dim, entities, tag, tags, remove, assertequaldim,
    boundary, skeleton, internal_boundary, exterior_boundary

end