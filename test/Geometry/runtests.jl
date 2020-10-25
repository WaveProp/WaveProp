using Test
using SafeTestsets

@safetestset "HyperRectangle tests" begin include("hyperrectangle_test.jl") end

@safetestset "Domain tests" begin include("domain_test.jl") end

@safetestset "Element tests" begin include("element_test.jl") end

@safetestset "Extruded element tests" begin include("extrudedElement_test.jl") end