using Test
using SafeTestsets

@safetestset "Geometry" begin include("Geometry/runtests.jl") end

@safetestset "Integration" begin include("Integration/runtests.jl") end

@safetestset "BIE" begin include("BIE/runtests.jl") end

@safetestset "IO" begin include("IO/runtests.jl") end

