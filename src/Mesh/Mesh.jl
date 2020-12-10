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
using WaveProp.Integration

import WaveProp.Geometry: ambient_dimension, geometric_dimension, domain, range

export GenericMesh, ElementIterator, nodes, etypes, elements, ent2tags, elements, near_interaction_list,
compute_quadrature!,SubMesh, AbstractMesh, CartesianMesh, meshgen

include("meshes.jl")
include("cartesianmesh.jl")
include("meshgen.jl")

end