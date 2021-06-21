using Test
using LinearAlgebra
using WaveProp
using WaveProp.Geometry
using WaveProp.Integration

@testset "Trapezoidal quadrature" begin
    q       = Trapezoidal{10}()
    D       = domain(q)
    @test D == ReferenceLine()
    x,w = q()
    @test sum(w) ≈ 1
    # integrate a periodic function. Should be very accurate.
    @test isapprox(integrate(x->cos(2π*x),q),0,atol=1e-10)
    @test integrate(x->sin(2π*x[1])^2,q) ≈ 0.5
end

@testset "TrapezoidalP quadrature" begin
    q = TrapezoidalP{10}()
    D = domain(q)
    @test D == ReferenceLine()
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
    d = ReferenceTriangle()
    # exact value for x^a*y^b integrate over reference triangle
    exa = (a,b) -> factorial(a)*factorial(b)/factorial(a+b+2)
    # check all quadrature implemented
    orders = keys(Integration.TRIANGLE_GAUSS_ORDER_TO_NPTS)
    for p in orders
        q     = Gauss(;domain=d,order=p)
        x,w   = q()
        for i in 0:p
            for j in 0:p-i
                @test integrate(x->x[1]^i*x[2]^j,q) ≈ exa(i,j)
            end
        end
    end
end

@testset "Gauss quad on tetrahedron" begin
    d = ReferenceTetrahedron()
    for p in (1,2)
        q = Gauss(;domain=d,order=p)
        x,w = q()
        @test sum(w) ≈ 1/6
    end
    # FIXME: check that we integrate all monomials up to `order` like in the
    # reference triangle
end

@testset "Tensor product quad on square" begin
    N,M = 10,12
    qx  = GaussLegendre(N)
    qy  = Fejer(M)
    q   = TensorProductQuadrature(qx,qy)
    a,b = 2*N-1,M-1 # maximum integration order of monomials
    @test integrate( x -> 1,q) ≈ 1
    f = x -> x[1]^a*x[2]^b
    @test integrate(f,q) ≈ 1/(a+1)*1/(b+1)
end
