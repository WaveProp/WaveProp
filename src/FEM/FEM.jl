module FEM

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using SparseArrays

# import all methods in WaveProp.INTERFACE_METHODS
using WaveProp
WaveProp.@import_interface

using WaveProp.Mesh
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation

# include("elementary_matrix.jl")
# include("dof_numbering.jl")

# export elementary_matrix
# export DOFNumbering, LocalDOFNumbering, dofs, elt2dof, local_dof_numbering
# export assembly

end
