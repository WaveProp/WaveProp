module IO

    using DocStringExtensions
    using GmshSDK
    using StaticArrays
    using Printf
    
    using WaveProp.Geometry    
    using WaveProp.Mesh
    using WaveProp.Utils

    export read_geo
    export read_msh
    
    include("gmshIO.jl")    
    
end