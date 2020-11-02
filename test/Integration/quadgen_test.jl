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
    x,w   = qrule(el)
    @test sum(w) ≈ 5
end

@testset "Triangle quadrature" begin
    qrule = Gauss{ReferenceTriangle,1}()
    vtx   = SVector(Point(0.,0.),Point(1.,0.),Point(0.,1.))
    F     = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
    x,w   = qrule(F)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    vtx = SVector(Point(-1.,0.),Point(1.,0.),Point(0.,1.))
    F   = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
    x,w = qrule(F)
    @test sum(w) ≈ 1 
end

@testset "Triangle surface quadrature" begin
    D     = ReferenceTriangle
    qrule = Gauss{D,1}()
    vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.))
    F   = LagrangeElement{D,3,3,Float64}(vtx)
    x,w = qrule(F)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    vtx = SVector(Point(-1.,0.,1.),Point(1.,0.,1.),Point(0.,1.,1.))
    F   = LagrangeElement{D,3,3,Float64}(vtx)
    x,w = qrule(F)
    @test sum(w) ≈ 1 
end

@testset "Tetrahedron quadrature" begin
    D     = ReferenceTetrahedron
    qrule = Gauss{D,1}()
    vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.),Point(0.,0.,1.))
    F   = LagrangeElement{D,4,3,Float64}(vtx)
    x,w = qrule(F)
    @test sum(w) ≈ 1/6
    # dilate by 2x and translate by 1 along  the tetrahedron
    vtx = SVector(Point(1.,0.,0.),Point(3.,0.,0.),Point(1.,2.,0.),Point(1.,0.,2.))
    F   = LagrangeElement{D,4,3,Float64}(vtx)
    x,w = qrule(F)
    @test sum(w) ≈ 1/6*2^3 
end
