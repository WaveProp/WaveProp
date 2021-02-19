module IO

using DocStringExtensions
using StaticArrays
using Printf
using RecipesBase
using WriteVTK
using OrderedCollections

using WaveProp.Utils    
using WaveProp.Geometry    
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.Nystrom

import WaveProp.Mesh: reference_nodes    

export read_geo
export read_msh
export vtk_mesh_file

include("gmshIO.jl")
include("plotsIO.jl")    
include("vtkIO.jl")

end # module
