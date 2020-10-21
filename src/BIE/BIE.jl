module BIE

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using IterativeSolvers

using WaveProp.Geometry
using WaveProp.Integration 
using WaveProp.Mesh

import WaveProp.Geometry: geometric_dimension, ambient_dimension, domain, range

export Laplace, Helmholtz, SingleLayerKernel, DoubleLayerKernel, IntegralOperator, NystromMesh, Density, γ₀, γ₁, SingleLayerPotential, DoubleLayerPotential, SingleLayerOperator, DoubleLayerOperator

include("pde.jl")
include("kernels.jl")
include("laplace.jl")
include("helmholtz.jl")
include("density.jl")
include("potential.jl")
include("integraloperators.jl")

# include("traits.jl")
# include("kernels.jl")
# include("helmholtz.jl")
# include("integraloperator.jl")

end # module