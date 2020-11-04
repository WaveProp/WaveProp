using Test
using LinearAlgebra
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.IO

@testset "Trapezoidal quadrature" begin
    q = Trapezoidal{10}()
    x,w = q()
    @test sum(w) ≈ 1
    # integrate a periodic function. Should be very accurate.
    @test isapprox(integrate(x->cos(2π*x[1]),q),0,atol=1e-10)
    @test integrate(x->sin(2π*x[1])^2,q) ≈ 0.5
end

@testset "Fejer quadrature" begin
    N = 5
    q = Fejer{N}()
    x,w = q()
    @test sum(w) ≈ 1
    # integrate all polynomial of degree N-1 exactly
    for n in 1:(N-1)
        @test integrate(x->x[1]^n,q) ≈ 1/(n+1)
    end
end

@testset "Gauss-Legendre quadrature" begin
    N = 5
    q = GaussLegendre{N}()
    x,w = q()
    @test sum(w) ≈ 1
    # integrate all polynomial of degree 2N-1 exactly
    for n in 1:(2*N-1)
        @test integrate(x->x[1]^n,q) ≈ 1/(n+1)
    end
end

@testset "Gauss quad on triangle" begin
    D= ReferenceTriangle
    for N in (1,3)
        q = Gauss{D,N}()
        x,w = q()
        @test sum(w) ≈ 1/2
    end
    # FIXME: check that we integrate all monomials up to `order`
end

@testset "Gauss quad on tetrahedron" begin
    D= ReferenceTetrahedron
    for N in (1,4)
        q = Gauss{D,N}()
        x,w = q()
        @test sum(w) ≈ 1/6
    end
    # FIXME: check that we integrate all monomials up to `order`
end

@testset "Tensor product quad on square" begin
    D   = ReferenceSquare
    N,M = 10,12
    qx  = GaussLegendre{N}()
    qy  = Fejer{M}()
    q   = TensorProductQuadrature(qx,qy)
    a,b = 2*N-1,M-1 # maximum integration order of monomials
    @test integrate(x->1,q) ≈ 1
    @test integrate(x->x[1]^a*x[2]^b,q) ≈ 1/(a+1)*1/(b+1)
end

@testset "Line quadrature" begin
    el     = LagrangeLine((1.,1.),(5.,4.))
    qrule  = GaussLegendre{1}()
    x,w   = qrule(el)
    @test sum(w) ≈ 5
end

@testset "Triangle quadrature" begin
    qrule = Gauss{ReferenceTriangle,1}()
    F     = LagrangeTriangle((0.,0.),(1.,0),(0.,1.))
    x,w   = qrule(F)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    F   = LagrangeTriangle((-1.,0),(1.,0),(0.,1.))
    x,w = qrule(F)
    @test sum(w) ≈ 1 
end

@testset "Triangle surface quadrature" begin
    qrule = Gauss{ReferenceTriangle,1}()
    F   = LagrangeTriangle((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
    x,w = qrule(F)
    @test sum(w) ≈ 1/2
    ## equilateral triangle
    F   = LagrangeTriangle((-1.,0.,1.),(1.,0.,1.),(0.,1.,1.))
    x,w = qrule(F)
    @test sum(w) ≈ 1 
end

@testset "Tetrahedron quadrature" begin
    D     = ReferenceTetrahedron
    qrule = Gauss{D,1}()
    F   = LagrangeTetrahedron((0,0,0.),(1.,0,0),(0,1,0),(0.,0.,1.))
    x,w = qrule(F)
    @test sum(w) ≈ 1/6
    # dilate by 2x and translate by 1 along  the tetrahedron
    F   = LagrangeTetrahedron((1,0,0),(3,0,0),(1,2,0),(1,0,2))
    x,w = qrule(F)
    @test sum(w) ≈ 1/6*2^3 
end
