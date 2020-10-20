"""
Julia library for solving PDEs with a focus on wave-propagation problems.

This package contains the following sub-modules:

- [`WaveProp.Geometry`](@ref)
- [`WaveProp.Integration`](@ref)
- [`WaveProp.FEM`](@ref)
- [`WaveProp.BIE`](@ref)

The exported names are:
$(EXPORTS)
"""
module WaveProp

using DocStringExtensions

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Integration/Integration.jl")

include("BIE/BIE.jl")

include("FEM/FEM.jl")

include("DDM/DDM.jl")

include("IO/IO.jl")

end # module
