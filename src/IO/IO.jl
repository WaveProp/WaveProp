module IO

    using DocStringExtensions
    using GmshSDK
    using StaticArrays
    
    using WaveProp.Geometry    
    using WaveProp.Utils

    export read_geo
    export read_msh
    
    include("gmshIO.jl")    
    
end