"""
    module Integration

Methods for integrating over instances of [`AbstractReferenceShape`](@ref).

Besides some standard quadrature rules used in the `FEM` module, also defines
singular integration routines useful for (weakly) singular integrands.

Exports:
$(EXPORTS)
"""
module Integration

using DocStringExtensions
using StaticArrays
using QuadGK
using LinearAlgebra

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Interpolation

# import all methods in WaveProp.INTERFACE_METHODS
using WaveProp
WaveProp.@import_interface

export
    # abstract types
    AbstractQuadratureRule,
    AbstractSingularityHandler,
    # types
    Gauss,
    Trapezoidal,
    TrapezoidalP,
    Fejer,
    GaussLegendre,
    TensorProductQuadrature,
    Kress,
    KressP,
    IMT,
    Duffy,
    SingularQuadratureRule,
    # functions
    integrate,
    qnodes,
    qweights,
    qnormals,
    singular_weights,
    singular_quadrature,
    TensorProductSingularityHandler

include("quadrulestables.jl")
include("quadrule.jl")
include("singularityhandler.jl")
include("singularquadrule.jl")

end
