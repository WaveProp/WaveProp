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
    q   = TensorProduct(qx,qy)
    a,b = 2*N-1,M-1 # maximum integration order of monomials
    @test integrate(x->1,q) ≈ 1
    @test integrate(x->x[1]^a*x[2]^b,q) ≈ 1/(a+1)*1/(b+1)
end

@testset "IMT on line" begin
    N      = 20
    qrule  = GaussLegendre{N}()
    x̂,ŵ    = qrule()
    e1 = integrate(log,x̂,ŵ) + 1 |> abs
    cov = IMT()
    x,w = push_forward_quad(cov,qrule)
    e2  = integrate(log,x,w) + 1 |> abs
    @test e2 < e1
end

@testset "Duffy on triangle" begin
    D   = ReferenceTriangle
    N   = 10
    cov = Duffy{2}()
    qx  = GaussLegendre{N}()
    qy     = qx
    qrule  = TensorProduct(qx,qy)
    x̂,ŵ = qrule()
    e1  = integrate(x̂,ŵ) do x
        d = sqrt((x[1]-1)^2 + x[2]^2)
        log(d)
    end       
    x,w   = push_forward_quad(cov,qrule)
    e2 = integrate(x,w) do x
        d = sqrt((x[1]-1)^2 + x[2]^2)
        log(d)
    end        
    # TODO: compare e1 and e2 to a reference value from e.g. quadgk. The value
    # of e2 should be much closer if the  change of variables was correct. 
end







