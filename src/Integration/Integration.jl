"""
    Methods for integrating [`AbstractElement`](@ref)
    
Besides some standard quadrature rules used in the [`FEM`](@ref) module, we also define singular integration routines for the [`BIE`](@ref) module.

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

import WaveProp.Geometry: domain, range, ambient_dimension, geometric_dimension

export integrate, Gauss, Trapezoidal, TrapezoidalP, Fejer, GaussLegendre, GenericQuadrature, TensorProductQuadrature, AbstractQuadratureRule, qnodes, qweights, qnormals

include("quadrature_rules.jl")

end
