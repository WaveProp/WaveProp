using Test
using SafeTestsets

@safetestset "Domain tests" begin include("domain_test.jl") end

@safetestset "Element tests" begin include("element_test.jl") end