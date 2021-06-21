using Test
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Interpolation
using WaveProp.Mesh
using LinearAlgebra
using StaticArrays

@testset "Lagrange elements" begin
    @testset "LagrangeLine" begin
        ## line in 3d
        vtx = ((0.,0.,0.),(1.,1.,1.))
        l   = LagrangeLine(vtx)
        @test domain(l) == ReferenceLine()
        @test length(return_type(l)) == 3
        @test dimension(domain(l)) == 1
        # makes no sense to ask for normal here. Make sure error is thrown.
        @test_throws AssertionError normal(l,SVector(0.5)) == SVector(-1.,1.)/√2
        ## line in 2d
        a = (0.,0.)
        b = (1.,1.)
        l   = LagrangeLine((a,b))
        @test domain(l) == ReferenceLine()
        @test length(return_type(l)) == 2
        @test dimension(domain(l)) == 1
        @test normal(l,SVector(0.5)) == SVector(1.,-1.)/√2
    end
    @testset "LagrangeTriangle" begin
        # triangle in 2d
        vtx = SVector(SVector(0.,0.),SVector(0.,1.),SVector(-1.,0))
        t   = LagrangeTriangle(vtx)
        @test length(return_type(t)) == 2
        @test dimension(domain(t)) == 2
        # triangle in 3d
        vtx = SVector(SVector(0.,0.,0.),SVector(0.,1.,0.),SVector(-1.,0,0.))
        t   = LagrangeTriangle(vtx)
        @test length(return_type(t)) == 3
        @test dimension(domain(t)) == 2
        @test normal(t,SVector(0.1,0.1)) == SVector(0,0,1.)
    end
    @testset "Tetrahedron" begin
        # TODO: add tests
    end
end
