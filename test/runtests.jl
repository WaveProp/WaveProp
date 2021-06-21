using Test
using SafeTestsets

@safetestset "Geometry" begin include("Geometry/runtests.jl") end

@safetestset "Interpolation" begin include("Interpolation/runtests.jl") end

@safetestset "Integration" begin include("Integration/runtests.jl") end

@safetestset "Mesh" begin include("Mesh/runtests.jl") end

# @safetestset "FEM" begin include("FEM/runtests.jl") end

# @safetestset "Nystrom" begin include("Nystrom/runtests.jl") end

# @safetestset "IO" begin include("IO/runtests.jl") end
