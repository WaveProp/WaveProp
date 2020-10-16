module BIE

using DocStringExtensions
using StaticArrays
using LinearAlgebra

using WaveProp.Geometry
using WaveProp.Integration    

export Helmholtz, SingleLayerKernel, DoubleLayerKernel, IntegralOperator

include("traits.jl")
include("pde.jl")
include("kernels.jl")
include("helmholtz.jl")
include("integraloperator.jl")

end # module