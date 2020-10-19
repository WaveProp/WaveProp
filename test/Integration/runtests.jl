using Test
using SafeTestsets

@safetestset "Quadrature rule tests" begin include("quadratureRule_test.jl") end

@safetestset "Quadgen tests" begin include("quadgen_test.jl") end
