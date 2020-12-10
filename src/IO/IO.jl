module IO

    using DocStringExtensions
    using GmshSDK
    using StaticArrays
    using Printf
    using RecipesBase
    using WriteVTK
    
    using WaveProp.Geometry    
    using WaveProp.Mesh
    using WaveProp.Utils
    using WaveProp.Nystrom

    import WaveProp.Geometry: reference_nodes    
    
    export read_geo
    export read_msh
    export vtk_mesh_file

    include("gmshIO.jl")
    include("plotsIO.jl")    
    include("vtkIO.jl")

end # module
