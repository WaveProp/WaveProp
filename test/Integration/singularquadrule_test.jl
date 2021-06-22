using Test
using WaveProp
using StaticArrays
using WaveProp.Geometry
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
    Ie    = -1
    for shand in [IMT(), Kress(),KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        Ia    = integrate(f,q)
        @test isapprox(Ia,Ie,rtol=1e-5)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,qstd)
        @test !isapprox(Istd,Ie,rtol=1e-3)
    end
    # test that KressP() can handle singularity at both endpoints, and Kress() cannot
    f     = (x) -> log(x) + log(1-x)
    Ie    = -2
    q     = SingularQuadratureRule(qstd,KressP())
    Ia    = integrate(f,q)
    @test isapprox(Ia,Ie,rtol=1e-5)
    q     = SingularQuadratureRule(qstd,Kress())
    Ia    = integrate(f,q)
    @test !isapprox(Ia,Ie,rtol=1e-5)
end

@testset "Duffy" begin
    q1d   = GaussLegendre(5)
    qstd  = TensorProductQuadrature(q1d,q1d)
    duffy = Duffy()
    qsin  = SingularQuadratureRule(qstd,duffy)
    x,w  = qsin()
    # regular kernel
    k    = x -> cos(x[1])*sin(x[2])
    Ie   = 1/2*(sin(1)-cos(1))
    Ia   = integrate(k,qsin)
    @test isapprox(Ie,Ia,rtol=1e-5)
    # singular kernel at right vertex
    s    = SVector(1,0)
    k    = x -> 1/norm(x-s)
    Ia   = integrate(k,qsin)
    Ie   = acosh(sqrt(2))
    @test isapprox(Ie,Ia,rtol=1e-5)
    # singular kernel at wrong vertex (should not be accurate, see test below)
    s    = SVector(0,1)
    k    = x -> 1/norm(x-s)
    Ie   = acosh(sqrt(2))
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
    k = (x) -> cos(x[1])
    Ie   = sin(1)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-6)
    # singular kernel at origin
    k  = (x) -> 1/norm(x)
    Ie = 1/2*log(17 + 12*sqrt(2)) # Wolfram alpha is your friend ;-)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-2)
end
