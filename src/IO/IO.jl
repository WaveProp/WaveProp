module IO

using DocStringExtensions
using StaticArrays
using Printf
using RecipesBase
using WriteVTK
using OrderedCollections

using WaveProp
using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using WaveProp.Nystrom

WaveProp.@import_interface

export read_geo
export read_msh
export vtk_mesh_file

include("gmshIO.jl")
include("plotsIO.jl")
include("vtkIO.jl")

end # module
