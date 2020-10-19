using Test
using StaticArrays
using LinearAlgebra
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.IO

@testset "Line quadrature" begin
    vtx   = SVector(Point(1.,1.),Point(5.,4.))
    el     = LagrangeElement{ReferenceLine,2,2,Float64}(vtx)
    qrule = GaussLegendre{1}()
    x,w   = quadgen(el,qrule)
    @test sum(w) ≈ 5
end

@testset "Triangle quadrature" begin
    qrule = Gauss{ReferenceTriangle,1}()
    vtx   = SVector(Point(0.,0.),Point(1.,0.),Point(0.,1.))
    F     = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
    x,w   = quadgen(F,qrule)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    vtx = SVector(Point(-1.,0.),Point(1.,0.),Point(0.,1.))
    F   = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1 
end

@testset "Triangle surface quadrature" begin
    D     = ReferenceTriangle
    qrule = Gauss{D,1}()
    vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.))
    F   = LagrangeElement{D,3,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    vtx = SVector(Point(-1.,0.,1.),Point(1.,0.,1.),Point(0.,1.,1.))
    F   = LagrangeElement{D,3,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1 
end

@testset "Tetrahedron quadrature" begin
    D     = ReferenceTetrahedron
    qrule = Gauss{D,1}()
    vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.),Point(0.,0.,1.))
    F   = LagrangeElement{D,4,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1/6
    # dilate by 2x and translate by 1 along  the tetrahedron
    vtx = SVector(Point(1.,0.,0.),Point(3.,0.,0.),Point(1.,2.,0.),Point(1.,0.,2.))
    F   = LagrangeElement{D,4,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1/6*2^3 
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        (lx,ly,lz) = widths = (1.,1.,2.)
        Ω, M  = WaveProp.IO.gmsh_box(;widths=widths)
        qrule = Gauss{ReferenceTriangle,1}()
        Q = quadgen(M,qrule;dim=2)
        A = 2*(lx*ly + lx*lz + ly*lz)
        @test A ≈ sum(Q.weights)
        qrule = Gauss{ReferenceTetrahedron,1}()
        Q = quadgen(M,qrule;dim=3)
        V = prod(widths)
        @test V ≈ sum(Q.weights)
    end
    @testset "Sphere" begin
        r = 0.5
        Ω, M = WaveProp.IO.gmsh_sphere(;radius=r)
        qrule = Gauss{ReferenceTriangle,1}()
        Q = quadgen(M,qrule;dim=2)
        A = 4π*r^2
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(Q.weights);atol=0.1)        
        qrule = Gauss{ReferenceTetrahedron,1}()
        Q = quadgen(M,qrule;dim=3)
        V = (4/3)*π*r^3
        @test isapprox(V,sum(Q.weights);atol=0.1)
    end
end