using WaveProp
using WaveProp.Geometry
using WaveProp.IO

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.geo")
    Ω = read_geo(fname)
    @test length(Ω) == 1
    @test geometric_dimension(Ω) == 3

    Γi = internal_boundary(Ω)
    Γe = external_boundary(Ω)

    @test length(Γi) == 0
    @test length(Γe) == 1
    @test geometric_dimension(Γe) == 2

    Ce = external_boundary(Γe)
    Ci = internal_boundary(Γe)
    S  = skeleton(Γe)
    @test length(Ce) == 2
    @test geometric_dimension(Ce) == 1
end

@testset "Sphere msh" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.msh")
    Ω, M = read_msh(fname)
    @test etypes(M) == [LagrangeLine{2,3,Float64},LagrangeTriangle{3,3,Float64},LagrangeTetrahedron{4,3,Float64},Point{3,Float64}] # mesh composed of gmsh simplices
    # Test internal creation of sphere
    Ω, M = WaveProp.IO.gmsh_sphere()    
    @test etypes(M) == [LagrangeLine{2,3,Float64},LagrangeTriangle{3,3,Float64},LagrangeTetrahedron{4,3,Float64},Point{3,Float64}] # mesh composed of gmsh simplices
end
