using Test
using SafeTestsets

@safetestset "Geometry" begin include("Geometry/runtests.jl") end

@testset "IO" begin include("Integration/runtests.jl") end

@testset "IO" begin include("IO/runtests.jl") end