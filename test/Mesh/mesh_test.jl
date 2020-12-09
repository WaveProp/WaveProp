using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Element iterator" begin
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    idx  = 2
    E    = etypes(M)[idx]
    iter = ElementIterator{E}(M)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Sub mesh" begin
    (lx,ly,lz) = widths = (1.,1.,2.)
    Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
    subM  = SubMesh(M,external_boundary(Ω))
    idx  = 2
    E    = etypes(M)[idx]
    iter = ElementIterator{E}(subM)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

