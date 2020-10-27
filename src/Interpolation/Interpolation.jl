module Interpolation
    
    using DocStringExtensions
    using StaticArrays

    export barycentric_lagrange_matrix, barycentric_lagrange_weights

    include("lagrange.jl")
    
end