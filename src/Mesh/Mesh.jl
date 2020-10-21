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

using WaveProp.Geometry
using WaveProp.Integration

import WaveProp.Geometry: ambient_dimension, geometric_dimension, domain, range
import WaveProp.Integration: quadgen

export GenericMesh, NystromMesh, etypes, ElementIterator

include("mesh.jl")

end