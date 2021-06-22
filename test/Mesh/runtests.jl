using SafeTestsets

@safetestset "Element tests" begin include("element_test.jl") end

@safetestset "CartesianMesh tests" begin include("cartesianmesh_test.jl") end
