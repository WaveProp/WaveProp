using Test
using StaticArrays
using LinearAlgebra
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.IO

# @testset "Line quadrature" begin
#     el  = ReferenceLine()
#     for N = 1:3
#         qrule = Integration.Gauss{N}()
#         x,w = quadgen(el,qrule)
#         @test sum(w) == 1
#     end
# end

# @testset "Triangle quadrature" begin
#     el  = ReferenceTriangle()
#     for N = 1:2
#         qrule = Integration.Gauss{N}()
#         x,w = quadgen(el,qrule)
#         @test sum(w) == 1/2
#     end
# end

# @testset "Tetrahedron quadrature" begin
#     el  = ReferenceTetrahedron()
#     for N = 1:2
#         qrule = Integration.Gauss{N}()
#         x,w = quadgen(el,qrule)
#         @test sum(w) == 1/6
#     end
# end

# @testset "Line quadrature" begin
#     vtx = SVector(Point(1.,1.),Point(5.,4.))
#     F   = LagrangeElement{ReferenceLine,2,2,Float64}(vtx)
#     qrule = Gauss{1}()
#     x,w = quadgen(F,qrule)
#     @test sum(w) ≈ 5
# end

# @testset "Triangle quadrature" begin
#     qrule = Gauss{1}()
#     vtx = SVector(Point(0.,0.),Point(1.,0.),Point(0.,1.))
#     F   = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
#     x,w = quadgen(F,qrule)
#     @test sum(w) ≈ 1/2
#     ## equilateral triangle
#     vtx = SVector(Point(-1.,0.),Point(1.,0.),Point(0.,1.))
#     F   = LagrangeElement{ReferenceTriangle,3,2,Float64}(vtx)
#     x,w = quadgen(F,qrule)
#     @test sum(w) ≈ 1 
# end

# @testset "Triangle surface quadrature" begin
#     qrule = Gauss{1}()
#     vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.))
#     F   = LagrangeElement{ReferenceTriangle,3,3,Float64}(vtx)
#     x,w = quadgen(F,qrule)
#     @test sum(w) ≈ 1/2
#     ## equilateral triangle
#     vtx = SVector(Point(-1.,0.,1.),Point(1.,0.,1.),Point(0.,1.,1.))
#     F   = LagrangeElement{ReferenceTriangle,3,3,Float64}(vtx)
#     x,w = quadgen(F,qrule)
#     @test sum(w) ≈ 1 
# end

@testset "Tetrahedron quadrature" begin
    qrule = Gauss{2}()
    vtx = SVector(Point(0.,0.,0.),Point(1.,0.,0.),Point(0.,1.,0.),Point(0.,0.,1.))
    F   = LagrangeElement{ReferenceTetrahedron,4,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1/6
    # dilate by 2x and translate by 1 along  the tetrahedron
    vtx = SVector(Point(1.,0.,0.),Point(3.,0.,0.),Point(1.,2.,0.),Point(1.,0.,2.))
    F   = LagrangeElement{ReferenceTetrahedron,4,3,Float64}(vtx)
    x,w = quadgen(F,qrule)
    @test sum(w) ≈ 1/6*2^3 
end

@testset "Mesh quadrature" begin
    @testset "Cube" begin
        fname = "/home/lfaria/projects/WaveProp/test/IO/cube.msh"
        Ω, M = read_msh(fname)
        qrule = Gauss{1}()
        Q = quadgen(M,qrule;dim=2)
        l = 1
        A = 6*l^2
        @test A ≈ sum(Q.weights)
        Q = quadgen(M,qrule;dim=3)
        V = l^3
        @test V ≈ sum(Q.weights)
    end
    @testset "Sphere" begin
        fname = "/home/lfaria/projects/WaveProp/test/IO/sphere_refined.msh"
        Ω, M = read_msh(fname)
        qrule = Gauss{1}()
        Q = quadgen(M,qrule;dim=2)
        r = 0.5
        A = 4π*r^2
        # the coarse tolerance below is because we use flat elements to
        # approximate the surface area and volume of a sphere
        @test isapprox(A,sum(Q.weights);atol=0.1)        
        Q = quadgen(M,qrule;dim=3)
        V = (4/3)*π*r^3
        @test isapprox(V,sum(Q.weights);atol=0.1)
    end
end

