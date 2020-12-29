module Nystrom

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Base.Threads
using SpecialFunctions

using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration 
using WaveProp.Integration 
using WaveProp.Mesh
using WaveProp.PhysicalProblem

import WaveProp.Geometry: geometric_dimension, ambient_dimension, domain, range
import WaveProp.Integration: singular_weights, qnodes, qweights, qnormals
import WaveProp.Mesh: etypes, elements

export Laplace, Helmholtz, SingleLayerKernel, DoubleLayerKernel, IntegralOperator, NystromMesh, Density, γ₀, γ₁, SingleLayerPotential, DoubleLayerPotential, SingleLayerOperator, DoubleLayerOperator, GreensCorrection, singular_weights_dim, single_double_layer, isinside, NystromMesh, dof

include("nystrommesh.jl")
include("greensformulae.jl")
include("lebedev.jl")
include("kernels.jl")
include("laplace.jl")
include("helmholtz.jl")
include("density.jl")
include("potential.jl")
include("integraloperators.jl")
include("dim.jl")
include("ldim.jl")
include("assemble.jl")
include("martensen_kussmaul.jl")

end # module
