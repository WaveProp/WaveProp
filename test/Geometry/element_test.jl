using Test
using WaveProp.Geometry
using StaticArrays

@testset "Basic shapes" begin
    @testset "Line" begin
        l = ReferenceLine()
        @test ambient_dimension(l) == 1
        @test geometric_dimension(l) == 1
        x = Point(0.5)
        @test x ∈ l
        x = Point(1.0)
        @test x ∈ l
        x = Point(1.1)
        @test !in(x,l) 
    end
    @testset "Triangle" begin
        t = ReferenceTriangle()
        @test ambient_dimension(t) == 2
        @test geometric_dimension(t) == 2
        x = Point(0.5,0.5)
        @test x ∈ t
        x = Point(1.0,0.0)
        @test x ∈ t
        x = Point(1.1,0.0)
        @test !in(x,t) 
    end
    @testset "Tetrahedron" begin
    t = ReferenceTetrahedron()
    @test ambient_dimension(t) == 3
    @test geometric_dimension(t) == 3
    x = Point(0.5,0.5,0.0)
    @test x ∈ t
    x = Point(1.0,0.0,0.0) # point on edge
    @test x ∈ t
    x = Point(1.1,0.0,0.0)
    @test !in(x,t) 
    end
end

@testset "Lagrange elements" begin
    @testset "LagrangeLine" begin
        ## line in 3d
        vtx = SVector(Point(0.,0.,0.),Point(1.,1.,1.))
        l   = LagrangeLine(vtx)
        @test domain(l) == ReferenceLine()
        @test ambient_dimension(l) == 3
        @test geometric_dimension(l) == 1
        # makes no sense to ask for normal here. Make sure error is thrown.
        @test_throws AssertionError normal(l,Point(0.5)) == Point(-1.,1.)/√2
        ## line in 2d
        vtx = SVector(Point(0.,0.),Point(1.,1.))
        l   = LagrangeLine(vtx)
        @test domain(l) == ReferenceLine()
        @test ambient_dimension(l) == 2
        @test geometric_dimension(l) == 1
        @test normal(l,Point(0.5)) == Point(-1.,1.)/√2
    end
    @testset "LagrangeTriangle" begin
        # triangle in 2d
        vtx = SVector(Point(0.,0.),Point(0.,1.),Point(-1.,0))
        t   = LagrangeTriangle(vtx)
        @test ambient_dimension(t) == 2
        @test geometric_dimension(t) == 2
        # triangle in 3d
        vtx = SVector(Point(0.,0.,0.),Point(0.,1.,0.),Point(-1.,0,0.))
        t   = LagrangeTriangle(vtx)
        @test ambient_dimension(t) == 3
        @test geometric_dimension(t) == 2
        @test normal(t,Point(0.1,0.1)) == Point(0,0,1.)
    end
    @testset "Tetrahedron" begin
        # TODO: add tests
    end
end