"""
    module PhysicalProblem
"""
module PhysicalProblem
    
    using DocStringExtensions

    using StaticArrays    

    import WaveProp.Geometry: ambient_dimension, geometric_dimension

    export 
        AbstractPDE, 
        Laplace, 
        Helmholtz, 
        Elastostatic, 
        default_density_eltype, 
        default_kernel_eltype

    include("pde.jl")    

end