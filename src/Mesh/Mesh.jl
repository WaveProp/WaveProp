"""
    module Mesh

Mesh data structures and the interface expected from an `AbstractMesh`.

Exports: $(EXPORTS)
"""

module Mesh 

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements
using OrderedCollections

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration

import WaveProp.Geometry: ambient_dimension, geometric_dimension, domain, range, boundary, parametrization
import WaveProp.Integration: qnodes, qweights, qnormals, integrate, normal, jacobian, singular_quadrature

export GenericMesh, ElementIterator, nodes, etypes, elements, ent2tags, el2qnodes, elements, near_interaction_list, nodes,
compute_quadrature!,SubMesh, AbstractMesh, CartesianMesh, meshgen, NystromMesh, dof, AbstractElement, LagrangeLine, LagrangeTriangle, LagrangeRectangle, LagrangeTetrahedron, LagrangeElement, ParametricElement, derivative, derivative2, measure, convert_to_2d, decompose

include("element.jl")
include("meshes.jl")
include("submesh.jl")
include("cartesianmesh.jl")
include("meshgen.jl")
include("decompose.jl")
include("elementiterator.jl")

end
