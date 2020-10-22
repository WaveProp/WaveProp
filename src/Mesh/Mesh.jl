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
import WaveProp.Integration: quadgen

export GenericMesh, etypes, ElementIterator, nodes, qnodes, qweights, qnormals, elements, nearest_element_list, compute_quadrature!

include("mesh.jl")
include("queries.jl")

end