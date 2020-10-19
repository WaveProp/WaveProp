"""
    Methods for integrating [`AbstractElement`](@ref)
    
Besides some standard quadrature rules used in the [`FEM`](@ref) module, we also define singular integration routines for the [`BIE`](@ref) module.

Exports:
$(EXPORTS)
"""
module Integration
    
using DocStringExtensions
using StaticArrays
using FastGaussQuadrature: gausslegendre
using LinearAlgebra

using WaveProp.Utils
using WaveProp.Geometry

export quadgen, Gauss, GenericQuadrature

include("quadratureRule.jl")
include("quadrature.jl")
include("quadgen.jl")

end