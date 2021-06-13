"""
    module Integration

Methods for integrating over instances of [`AbstractReferenceShape`](@ref).

Besides some standard quadrature rules used in the `FEM` module, this module also
defines singular integration routines useful for (weakly) singular integrands.

Exports:
$(EXPORTS)
"""
module Integration

using DocStringExtensions
using StaticArrays
using QuadGK
using LinearAlgebra

import QuadGK: quadgk
import HCubature: hcubature

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Interpolation

import WaveProp.Geometry:
    domain,
    range,
    ambient_dimension,
    geometric_dimension,
    jacobian

import WaveProp.Interpolation: barycentric_lagrange_weights

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

include("quadrule.jl")
include("quadgk.jl")
include("hcubature.jl")
include("singularityhandler.jl")
include("singularquadrule.jl")

end
