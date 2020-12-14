using SafeTestsets

@safetestset "Element tests" begin include("element_test.jl") end

@safetestset "Mesh tests" begin include("mesh_test.jl") end

@safetestset "CartesianMesh tests" begin include("cartesianmesh_test.jl") end

@safetestset "Meshgen tests" begin include("meshgen_test.jl") end