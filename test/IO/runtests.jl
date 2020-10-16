using Test
using WaveProp.Geometry
using WaveProp.IO

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.geo")
    Ω = read_geo(fname)
    @test length(Ω) == 1
    @test dim(Ω) == 3

    Γi = internal_boundary(Ω)
    Γe = external_boundary(Ω)

    @test length(Γi) == 0
    @test length(Γe) == 1
    @test dim(Γe) == 2

    Ce = external_boundary(Γe)
    Ci = internal_boundary(Γe)
    S  = skeleton(Γe)
    @test length(Ce) == 2
    @test dim(Ce) == 1
end

@testset "Sphere msh" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.msh")
    Ω, M = read_msh(fname)
    @test etypes(M) == [LagrangeLine{2},LagrangeTriangle{3},LagrangeTetrahedron{4},Point{3,Float64}] # mesh composed of gmsh simplices
    # Test internal creation of sphere
    Ω, M = WaveProp.IO.gmsh_sphere()    
    @test etypes(M) == [LagrangeLine{2},LagrangeTriangle{3},LagrangeTetrahedron{4},Point{3,Float64}] # mesh composed of gmsh simplices
end

