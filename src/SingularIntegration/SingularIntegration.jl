"""
    module SingularIntegration

Implement specialised integration routines for non-smooth integrands. In
particular, this module focuses on integrands which typically arise in boundary
integral equations.
"""
module SingularIntegration

using DocStringExtensions
using StaticArrays
using QuadGK

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation

import WaveProp.Geometry: domain, range, jacobian
import WaveProp.Integration: integrate

import QuadGK: quadgk

export Kress, KressP, IMT, Window, SingularQuadratureRule, singular_weights, singular_quadrature, Duffy, TensorProductSingularityHandler

include("singularityhandler.jl")
include("singularquadrule.jl")

end # module
