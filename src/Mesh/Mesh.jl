"""
    Mesh

This module defines various meshes used, as well as the interface expected from
an `AbstractMesh`.

Exports:
$(EXPORTS)
"""

module Mesh 

using DocStringExtensions
using StaticArrays
using LinearAlgebra

using WaveProp.Utils
using WaveProp.Geometry

import WaveProp.Geometry: ambient_dimension, geometric_dimension, domain, range

export GenericMesh, ElementIterator, nodes, etypes, el2nodes, ent2tags, elements, near_interaction_list,
compute_quadrature!,SubMesh, AbstractMesh

include("meshes.jl")
include("queries.jl")

end