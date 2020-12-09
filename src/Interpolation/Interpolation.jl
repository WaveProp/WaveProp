module Interpolation
    
    using DocStringExtensions
    using StaticArrays
    using GmshSDK

    using WaveProp.Utils    
    using WaveProp.Geometry    

    export barycentric_lagrange_matrix, barycentric_lagrange_weights

    include("lagrange.jl")
    
end