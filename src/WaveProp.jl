"""
Julia library for solving PDEs with a focus on wave-propagation problems.
"""
module WaveProp

using DocStringExtensions

include("interface.jl")

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Interpolation/Interpolation.jl")

# include("Integration/Integration.jl")

# include("Mesh/Mesh.jl")

# include("PhysicalProblem/PhysicalProblem.jl")

# include("Nystrom/Nystrom.jl")

# include("FEM/FEM.jl")

# include("IO/IO.jl")

# export all methods in INTERFACE_METHODS
@export_interface

end # module
