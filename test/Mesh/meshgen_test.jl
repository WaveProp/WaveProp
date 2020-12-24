using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Circle" begin
    geo = Circle()
    bnd = boundary(geo)
    Ω,mesh =  Mesh.meshgen(geo;gridsize=0.1)
    @test ambient_dimension(mesh) == 2
end