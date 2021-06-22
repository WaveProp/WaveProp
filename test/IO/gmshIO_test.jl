using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.IO
using WaveProp.Mesh

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.geo")
    Geometry.clear_entities!()
    Ω,M = read_geo(fname)
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
    Geometry.clear_entities!()
    fname = joinpath(@__DIR__,"sphere.msh")
    Ω, M  = read_msh(fname)
    etypes(M)
    T = SVector{3,Float64}
    @test etypes(M) == [LagrangeLine{2,T},LagrangeTriangle{3,T},LagrangeTetrahedron{4,T},T] # mesh composed of gmsh simplices
    # Test internal creation of sphere
    Ω, M = WaveProp.IO.gmsh_sphere()
    @test etypes(M) == [LagrangeLine{2,T},LagrangeTriangle{3,T},LagrangeTetrahedron{4,T},T] # mesh composed of gmsh simplices
end

@testset "Disk" begin
    # Test internal creation of disk
    Geometry.clear_entities!()
    Ω, M = WaveProp.IO.gmsh_disk()
    T = SVector{2,Float64}
    @test etypes(M) == [LagrangeLine{2,T},LagrangeTriangle{3,T},T] # mesh composed of gmsh simplices
end

@testset "Element iterator" begin
    Geometry.clear_entities!()
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    idx  = 2
    E    = etypes(M)[idx]
    iter = ElementIterator{E}(M)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Sub mesh" begin
    Geometry.clear_entities!()
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    subM  = SubMesh(M,external_boundary(Ω))
    idx  = 2
    E    = etypes(M)[idx]
    iter = ElementIterator{E}(subM)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end
