using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.IO
using WaveProp.Interpolation
using WaveProp.Mesh
using StaticArrays

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    fname = joinpath(@__DIR__,"sphere.geo")
    clear_entities!()
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
    clear_entities!()
    fname = joinpath(@__DIR__,"sphere.msh")
    Ω, M  = read_msh(fname)
    # Test internal creation of sphere
    Ω, M = WaveProp.IO.gmsh_sphere(dim=2)
end

@testset "Disk" begin
    # Test internal creation of disk
    clear_entities!()
    Ω, M = WaveProp.IO.gmsh_disk()
end

@testset "Element iterator" begin
    clear_entities!()
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    E    = keys(M) |> first
    iter = ElementIterator{E}(M)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Sub mesh" begin
    clear_entities!()
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    subM  = SubMesh(M,external_boundary(Ω))
    E    = keys(subM) |> first
    iter = ElementIterator(subM,E)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end
