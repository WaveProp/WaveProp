module FEM

using DocStringExtensions
using StaticArrays
using LinearAlgebra

using WaveProp
using WaveProp.Mesh
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation

include("elementary_matrix.jl")

export elementary_matrix

end