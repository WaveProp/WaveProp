module Interpolation
    
    using DocStringExtensions
    using StaticArrays
    using LinearAlgebra

    using WaveProp.Utils    
    using WaveProp.Geometry    

    export 
        barycentric_lagrange_matrix, 
        barycentric_lagrange_weights,
        MonomialBasis,
        vandermond,
        lagrange_basis,
        PolynomialBasis,
        LagrangeBasis, gradLagrangeBasis,
        trace, grad

    include("lagrange.jl")
    include("monomial.jl")
    include("polynomialbasis.jl")


    
end