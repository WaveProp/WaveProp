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
        vtx = SVector(Point(0.,0.,0.),Point(1.,1.,1.))
        l   = LagrangeLine{2}(vtx)
        @test reference_element(l) == ReferenceLine()
        @test ambient_dimension(l) == 3
        @test geometric_dimension(l) == 1
    end
    # @testset "LagrangeTriangle" begin
    #     t = ReferenceTriangle()
    #     @test ambient_dimension(t) == 2
    #     x = Point(0.5,0.5)
    #     @test x ∈ t
    #     x = Point(1.0,0.0)
    #     @test x ∈ t
    #     x = Point(1.1,0.0)
    #     @test !in(x,t) 
    # end
    # @testset "Tetrahedron" begin
    # t = ReferenceTetrahedron()
    # @test ambient_dimension(t) == 3
    # x = Point(0.5,0.5,0.0)
    # @test x ∈ t
    # x = Point(1.0,0.0,0.0) # point on edge
    # @test x ∈ t
    # x = Point(1.1,0.0,0.0)
    # @test !in(x,t) 
    # end
end