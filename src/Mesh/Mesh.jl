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
using WaveProp.Interpolation
using WaveProp.Integration

# import all methods in WaveProp.INTERFACE_METHODS
using WaveProp
WaveProp.@import_interface

export
    # abstract types
    AbstractMesh,
    # structs
    GenericMesh,
    ElementIterator,
    NodeIterator,
    SubMesh,
    UniformCartesianMesh,
    ParametricElement,
    # methods
    nodes,
    elements,
    ent2tags,
    el2qnodes,
    elements,
    near_interaction_list,
    nodes,
    compute_quadrature!,
    dof,
    derivative,
    derivative2,
    measure,
    convert_to_2d,
    decompose,
    mesh

include("abstractmesh.jl")
include("genericmesh.jl")
include("cartesianmesh.jl")
include("submesh.jl")
include("decompose.jl")

end
