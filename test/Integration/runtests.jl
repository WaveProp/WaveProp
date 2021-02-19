using Test
using SafeTestsets

@safetestset "Quadrature rules" begin include("quadrule_test.jl") end

@safetestset "Singular quadrature rules" begin include("singularquadrule_test.jl") end
