module BIE

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Base.Threads

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration 
using WaveProp.SingularIntegration 
using WaveProp.Mesh

import WaveProp.Geometry: geometric_dimension, ambient_dimension, domain, range
import WaveProp.SingularIntegration: singular_weights

export Laplace, Helmholtz, SingleLayerKernel, DoubleLayerKernel, IntegralOperator, NystromMesh, Density, γ₀, γ₁, SingleLayerPotential, DoubleLayerPotential, SingleLayerOperator, DoubleLayerOperator, GreensCorrection, singular_weights

include("pde.jl")
include("lebedev.jl")
include("kernels.jl")
include("laplace.jl")
include("helmholtz.jl")
include("density.jl")
include("potential.jl")
include("integraloperators.jl")
include("densityinterpolation.jl")
include("localDIM.jl")
include("singularweights.jl")

end # module
