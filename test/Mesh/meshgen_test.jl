using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Circle" begin
    geo = Circle()
    bnd = boundary(geo)
    @test length(bnd) == 1
    Ω    = Domain([geo])
    mesh =  Mesh.meshgen(Ω;h=0.1)
    @test ambient_dimension(mesh) == 2
end