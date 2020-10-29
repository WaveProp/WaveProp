module IO

    using DocStringExtensions
    using GmshSDK
    using StaticArrays
    using Printf
    using WriteVTK
    
    using WaveProp.Geometry    
    using WaveProp.Mesh
    using WaveProp.Utils
    
    include("gmshIO.jl")    
    include("vtkIO.jl")    

    export read_geo
    export read_msh
    
    export vtk_mesh_file
    
end