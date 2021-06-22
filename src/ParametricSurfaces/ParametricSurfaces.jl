module ParametricSurfaces

using DocStringExtensions
using StaticArrays
using LinearAlgebra
using ForwardDiff # for computing derivatives of parametric elements

using WaveProp
using WaveProp.Utils
using WaveProp.Geometry
using WaveProp.Interpolation
using WaveProp.Mesh

WaveProp.@import_interface

export
    #types
    ParametricEntity,
    ParametricElement,
    #functions
    line,
    meshgen

include("parametricentity.jl")
include("parametricelement.jl")
include("meshgen.jl")
include("simpleshapes.jl")

end # module
