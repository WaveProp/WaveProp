"""
    module SingularIntegration

Implement specialised integration routines for non-smooth integrands. In
particular, this module focuses on integrands which typically arise in boundary
integral equations.

Exports:
$(EXPORTS)
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

export Kress, IMT, SingularQuadratureRule, singular_weights, Duffy, TensorProductHandler

include("singularityhandler.jl")
include("singularquadrule.jl")

end # module
