using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh

@testset "Element iterator 1d" begin
    d = ReferenceLine()
    h = 0.1
    mesh = CartesianMesh(d;h=h)
    els = ElementIterator(mesh)
    @test length(els) == 10
    @test eltype(els) == HyperRectangle{1,Float64}
    @test first(els) == HyperRectangle((0,),(h,))
    for (n,el) in enumerate(els)
        @test el ≈ HyperRectangle((n-1)*h,n*h)
    end
end

@testset "Element iterator 2d" begin
    # d    = ReferenceSquare()
    # h    = 0.1
    # mesh = CartesianMesh(d;h=(h,h),n=(1,1))
    # els = ElementIterator(mesh)
    @test_broken length(els) == 100
    @test_broken eltype(els) == HyperRectangle{2,Float64}
    @test_broken first(els) == HyperRectangle((0,0),(h,h))
end