"""
    module PhysicalProblem
"""
module PhysicalProblem

using DocStringExtensions

using StaticArrays

# import all methods in WaveProp.INTERFACE_METHODS
using WaveProp
WaveProp.@import_interface

export
    AbstractPDE,
    Laplace,
    Helmholtz,
    Elastostatic,
    Maxwell,
    MaxwellCFIE,
    default_density_eltype,
    default_kernel_eltype,
    parameters

include("pde.jl")

end
