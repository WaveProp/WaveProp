module Nystrom

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Base.Threads
using SpecialFunctions
using OrderedCollections
using QuadGK

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
    MaxwellCFIE,
    SingleLayerKernel,
    DoubleLayerKernel,
    AdjointDoubleLayerKernel,
    HyperSingularKernel,
    GenericKernel,
    IntegralOperator,
    NystromMesh,
    Density,
    TangentialDensity,
    SingleLayerPotential,
    DoubleLayerPotential,
    SingleLayerOperator,
    DoubleLayerOperator,
    AdjointDoubleLayerOperator,
    HyperSingularOperator,
    # methods
    ncross,
    trace,
    coords,
    weight,
    γ₀,
    γ₁,
    singular_weights_dim,
    single_double_layer,
    isinside,
    dom2dof,
    assemble,
    qweights,
    diagonal_ncross_matrix

include("nystrommesh.jl")
include("abstractkernel.jl")
include("greensformulae.jl")
include("lebedev.jl")
include("laplace.jl")
include("helmholtz.jl")
include("elastostatic.jl")
include("maxwell.jl")
include("maxwellCFIE.jl")
include("density.jl")
include("potential.jl")
include("integraloperators.jl")
include("dim.jl")
include("maxwellCFIE_dim.jl")
include("assemble_gk.jl")
# include("ldim.jl")
# include("assemble.jl")
# include("martensen_kussmaul.jl")

end # module
