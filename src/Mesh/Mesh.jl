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

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration

import WaveProp.Geometry: ambient_dimension, geometric_dimension, domain, range
import WaveProp.Integration: quadgen, qnodes, qweights, qnormals

export GenericMesh, etypes, ElementIterator, nodes, elements, nearest_element_list, compute_quadrature!, el2qnodes

include("meshes.jl")
include("queries.jl")

end