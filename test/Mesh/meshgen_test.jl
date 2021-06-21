using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Disk" begin
    disk = Geometry.Circle() # abstract entity
    Γ = boundary(disk) |> Domain
    M = meshgen(Γ,(10,))
    # plot(M,Γ)
    @test entities(Γ) == boundary(disk)
    @test geometric_dimension(disk) == 2
end

@testset "Ball" begin
    ball = Geometry.Sphere() # abstract entity
    Γ = boundary(ball) |> Domain
    M = meshgen(Γ,(1,1))
    # plot(M,Γ)
    @test entities(Γ) == boundary(ball)
    @test geometric_dimension(ball) == 3
end
