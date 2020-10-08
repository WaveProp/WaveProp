"""
    Geometry

This module defines the basic geometrical objects manipulated. 

Exports:
$(EXPORTS)
"""
module Geometry 

import Base.==, Base.length, Base.in, Base.iterate
import Base.getindex, Base.lastindex, Base.isempty
import Base.union, Base.setdiff, Base.intersect, Base.issubset

using DocStringExtensions
using StaticArrays

include("point.jl")
include("domain.jl")
export ElementaryEntity, Domain,
    ==, length, in, iterate,
    isempty, union, setdiff, intersect, issubset, remove,
    dim, entities, tag, tags,
    assertequaldim, boundary

end