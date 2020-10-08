using Test
using SafeTestsets

@time @safetestset "Geometry" begin include("Geometry/runtests.jl") end