using Test
using WaveProp.Geometry
using WaveProp.Interpolation
using WaveProp.Integration
using LinearAlgebra

@testset "Reference segment" begin
    qstd = GaussLegendre(20)    
    d     = domain(qstd)    
    # simple test on smooth integrand
    for shand in [IMT(), Kress(), KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        @test isapprox(integrate(cos,q),sin(1),rtol=1e-2)
    end
    # non-smooth integrand, singularity at 0
    f     = (x) -> log(abs(x))
    I     = integrate(f,d;rtol=1e-16)
    for shand in [IMT(), Kress(),KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        Ia    = integrate(f,q)
        @test isapprox(Ia,I,rtol=1e-3)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,qstd)
        @test !isapprox(Istd,I,rtol=1e-3)
    end
    # non-smooth integrand, singularity somewhere inside
    s     = rand()
    f     = (x) -> x[1]==s ? 0.0 : log(abs(x[1]-s))
    I     = integrate(f,d;rtol=1e-8)
    for shand in [IMT(), Kress(),KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        x,w   = singular_quadrature(q,s)
        Ia    = integrate(f,x,w)
        @debug I - Ia
        @test isapprox(Ia,I,rtol=1e-3)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,q)
        @test !isapprox(Istd,I,rtol=1e-3)
    end
    # quadrature for integrating a factored kernel singular at location s
    k = (x) -> x[1]==s ? 0. : log(abs(x[1]-s))
    ϕ = (x) -> cos(x[1])
    f = (x) -> k(x)*ϕ(x)
    I = integrate(f,d)
    for shand in [IMT(), Kress(),KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        x,w   = singular_quadrature(k,q,s)
        Ia    = integrate(ϕ,x,w)
        @debug I - Ia
        @test isapprox(Ia,I,rtol=1e-5)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,q)
        @test !isapprox(Istd,I,rtol=1e-5)
    end
    # weights for integrating a factored kernel singular at location s
    x  = qnodes(qstd)
    for shand in [IMT(), Kress(),KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        ws    = singular_weights(k,qstd,q,s)
        Ia    = integrate(ϕ,x,ws)
        @test isapprox(Ia,I,rtol=1e-3)
    end
end

@testset "Duffy" begin
    q1d   = GaussLegendre(5)
    qstd  = TensorProductQuadrature(q1d,q1d)
    duffy = Duffy()
    qsin  = SingularQuadratureRule(qstd,duffy)
    x,w  = qsin()
    # regular kernel
    k    = x -> cos(x[1])*sin(x[2])
    Ie    = integrate(k,ReferenceTriangle)
    Ia   = integrate(k,qsin)
    @test isapprox(Ie,Ia,rtol=1e-5)
    # singular kernel at right vertex
    s    = SVector(1,0) 
    k    = x -> 1/norm(x-s)
    Ia   = integrate(k,qsin) 
    Ie    = integrate(k,ReferenceTriangle) # uses quadgk
    @test isapprox(Ie,Ia,rtol=1e-5)
    # singular kernel at wrong vertex (should not be accurate, see test below)
    s    = SVector(0,1) 
    k    = x -> 1/norm(x-s)
    Ie    = integrate(k,ReferenceTriangle)
    Ia   = integrate(k,qsin)
    @test !isapprox(Ie,Ia,rtol=1e-5)
end

@testset "2d Kress" begin
    q1d   = GaussLegendre{10}()
    qstd  = TensorProductQuadrature(q1d,q1d)
    s1d   = Kress(order=2)
    sing_handler = TensorProductSingularityHandler(s1d,s1d)
    qsin = SingularQuadratureRule(qstd,sing_handler)
    x,w  = qsin()    
    # regular kernel
    k(x) = cos(x[1])
    Ie    = integrate(k,ReferenceSquare)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-2)
    # singular kernel at origin
    k(x) = 1/norm(x)
    Ie    = integrate(k,ReferenceSquare)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-2)
end

