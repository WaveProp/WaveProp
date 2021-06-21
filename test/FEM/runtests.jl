using Test
using SafeTestsets

@safetestset "FEM Matrix assembly" begin include("fem_assembly.jl") end
