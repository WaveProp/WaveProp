using Test
using WaveProp.Geometry
using WaveProp.Integration
using LinearAlgebra

@testset "Reference segment" begin
    qstd = GaussLegendre(20)    
    # simple test on smooth integrand
    for shand in [IMT(), Kress(), KressP()]
        q     = SingularQuadratureRule(qstd,shand)
        @test isapprox(integrate(cos,q),sin(1),rtol=1e-2)
    end
    # non-smooth integrand
    for shand in [IMT(), Kress(),KressP()]
        d     = domain(qstd)    
        f     = (x) -> log(abs(x))
        q     = SingularQuadratureRule(qstd,shand)
        I     = integrate(f,d;rtol=1e-16)
        Ia    = integrate(f,q)
        @test isapprox(Ia,I,rtol=1e-3)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,qstd)
        @test !isapprox(Istd,I,rtol=1e-3)
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

# @testset "2d Kress" begin
#     q1d   = GaussLegendre{10}()
#     qstd  = TensorProductQuadrature(q1d,q1d)
#     s1d   = Kress(order=2)
#     sing_handler = TensorProductQuadratureHandler(s1d,s1d)
#     qsin = SingularQuadratureRule(qstd,sing_handler)
#     x,w  = qsin()    
#     # regular kernel
#     k(x) = cos(x[1])
#     I    = quadgk(k,ReferenceSquare)
#     Isin = sum(k.(x).*w)
#     @test isapprox(I,Isin,rtol=1e-2)
#     # singular kernel at origin
#     k(x) = 1/norm(x)
#     I    = quadgk(k,ReferenceSquare)
#     Isin = sum(k.(x).*w)
#     @test isapprox(I,Isin,rtol=1e-2)
#     # singular kernel elsewhere
#     s   = SVector(0.47,0.38)
#     x,w  = qsin(ReferenceSquare(),s)
#     k(x) = 1/norm(x-s)
#     I    = quadgk(k,ReferenceSquare)
#     Isin = sum(k.(x).*w)
#     @test isapprox(I,Isin,rtol=1e-2)
# end

