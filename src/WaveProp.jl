"""
Julia library for solving PDEs with a focus on wave-propagation problems.
"""
module WaveProp

using DocStringExtensions

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Mesh/Mesh.jl")

include("Interpolation/Interpolation.jl")

include("Integration/Integration.jl")

include("SingularIntegration/SingularIntegration.jl")

include("PhysicalProblem/PhysicalProblem.jl")

include("Nystrom/Nystrom.jl")

# include("FEM/FEM.jl")

# include("DDM/DDM.jl")

include("IO/IO.jl")

end # module
