module Nystrom

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Base.Threads
using SpecialFunctions
using OrderedCollections

using WaveProp
using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation
using WaveProp.Mesh
using WaveProp.PhysicalProblem

WaveProp.@import_interface

export
    # abstract types
    # types
    NystromMesh,
    NystromDOF,
    Laplace,
    Helmholtz,
    Elastostatic,
    Maxwell,
    Stokes,
    Maxwell,
    SingleLayerKernel,
    DoubleLayerKernel,
    AdjointDoubleLayerKernel,
    HyperSingularKernel,
    GenericKernel,
    IntegralOperator,
    NystromMesh,
    Density,
    SingleLayerPotential,
    DoubleLayerPotential,
    SingleLayerOperator,
    DoubleLayerOperator,
    AdjointDoubleLayerOperator,
    HyperSingularOperator,
    # methods
    coords,
    weight,
    γ₀,
    γ₁,
    singular_weights_dim,
    single_double_layer,
    isinside,
    dom2dof,
    assemble,
    qweights

include("nystrommesh.jl")
include("abstractkernel.jl")
include("greensformulae.jl")
include("lebedev.jl")
include("laplace.jl")
include("helmholtz.jl")
include("elastostatic.jl")
include("maxwell.jl")
include("density.jl")
include("potential.jl")
include("integraloperators.jl")
include("dim.jl")
# include("ldim.jl")
# include("assemble.jl")
# include("martensen_kussmaul.jl")

end # module
