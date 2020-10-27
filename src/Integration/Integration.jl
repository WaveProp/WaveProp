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

import QuadGK: quadgk

using WaveProp.Utils
using WaveProp.Geometry

import WaveProp.Geometry: domain, range, ambient_dimension, geometric_dimension

export integrate, quadgen, Gauss, Trapezoidal, Fejer, GaussLegendre, GenericQuadrature, TensorProduct, lebedev_points, qnodes, qweights, qnormals, push_forward_quad, push_forward_quad_with_normal

include("quadratureRule.jl")
include("lebedev.jl")
include("quadrature.jl")
include("quadgen.jl")
include("quadgk.jl")

end