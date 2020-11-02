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
import WaveProp.Integration: qnodes, qweights, qnormals

export GenericMesh, ElementIterator, nodes, etypes, el2nodes, ent2tags,
qnodes, qweights, qnormals, elements, nearest_element_list,
compute_quadrature!

include("meshes.jl")
include("queries.jl")

end