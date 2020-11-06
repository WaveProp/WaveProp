using Test
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.SingularIntegration
using LinearAlgebra

@testset "Reference segment" begin
    qstd = GaussLegendre(20)    
    # simple test on smooth integrand
    for shand in [IMT(), Kress()]
        q     = SingularQuadratureRule(qstd,shand)
        @test isapprox(integrate(cos,q),sin(1),rtol=1e-2)
    end
    # non-smooth integrand
    for shand in [IMT(), Kress()]
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

@testset "Curved lines" begin
    qstd = GaussLegendre(20)    
    el   = LagrangeLine((0,0),(1,1),(1/2,1/4))
    # simple test on smooth integrand
    f     = x->cos(x[1])
    I     = integrate(f,el)
    for shand in [IMT(), Kress(),Window()]
        q     = SingularQuadratureRule(qstd,shand)
        Ia    = integrate(f,q,el)
        @debug Ia - I
        @test isapprox(Ia,I,rtol=1e-6)
    end
    # non-smooth integrand
    f     = x->log(abs(x[1])*cos(x[1]))
    I     = integrate(f,el)
    for shand in [IMT(), Kress(),Window()]
        q     = SingularQuadratureRule(qstd,shand)    
        Ia    = integrate(f,q,el)
        @debug Ia - I
        @test isapprox(Ia,I,rtol=1e-5)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,qstd,el)
        @debug Istd - I
        @test !isapprox(Istd,I,rtol=1e-5)
    end
end

@testset "Curved line with singularity inside" begin
    qstd = GaussLegendre(20)    
    el   = LagrangeLine((0,0),(1,1),(1/2,1/4))
    vs   = 1/3
    xs   = el(vs)
    f    = x -> x==xs ? 0.0 : log(norm(x-xs)*cos(x[1]))
    I    = integrate(f,el)
    for shand in [IMT(), Kress(),Window()]
        q     = SingularQuadratureRule(qstd,shand)    
        x,w   = q(el,vs)
        Ia    = integrate(f,x,w)
        @debug I-Ia
        @test isapprox(I,Ia;rtol=1e-5)
        Ia    = integrate(f,q,el,vs)
        @test isapprox(I,Ia;rtol=1e-5)
        # check that the `naive` integration woudl have failed the test
        Istd    = integrate(f,qstd,el)
        @test !isapprox(Istd,I,rtol=1e-5)
    end
end

# @testset "Duffy" begin
#     @testset "Singularity at vertex" begin
#         q1d   = GaussLegendre(5)
#         qstd  = TensorProductQuadrature(q1d,q1d)
#         duffy = Duffy{2}()
#         qsin = SingularQuadratureRule(qstd,duffy)
#         x,w  = qsin()
#         # regular kernel
#         k(x) = cos(x[1])
#         I    = quadgk(k,ReferenceTriangle)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-5)
#         # singular kernel at right vertex
#         s    = Point(1,0) 
#         k(x) = 1/norm(x-s)
#         I    = quadgk(k,ReferenceTriangle)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-5)
#         # singular kernel at wrong vertex (should not be accurate, see test below)
#         s    = Point(0,1) 
#         k(x) = 1/norm(x-s)
#         I    = quadgk(k,ReferenceTriangle)
#         Isin = sum(k.(x).*w)
#         @test !isapprox(I,Isin,rtol=1e-5)
#     end

#     @testset "Singularity inside" begin
#         q1d   = GaussLegendre{5}()
#         qstd  = TensorProductQuadrature(q1d,q1d)
#         duffy = Duffy{2}()
#         qsin = SingularQuadratureRule(qstd,duffy)
#         s    = Point(0.37,0.25) 
#         # triangle case
#         x,w  = qsin(ReferenceTriangle(),s)
#         # regular kernel
#         k(x) = cos(x[1])
#         I    = quadgk(k,ReferenceTriangle)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-5)
#         # singular kernel at interior point
#         k(x) = 1/norm(x-s)
#         I    = quadgk(k,ReferenceTriangle)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-2)
#         # square case
#         x,w  = qsin(ReferenceSquare(),s)
#         # regular kernel
#         k(x) = cos(x[1])
#         I    = quadgk(k,ReferenceSquare)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-5)
#         # singular kernel at interior point
#         k(x) = 1/norm(x-s)
#         I    = quadgk(k,ReferenceSquare)
#         Isin = sum(k.(x).*w)
#         @test isapprox(I,Isin,rtol=1e-2)
#     end    
# end

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
#     s   = Point(0.47,0.38)
#     x,w  = qsin(ReferenceSquare(),s)
#     k(x) = 1/norm(x-s)
#     I    = quadgk(k,ReferenceSquare)
#     Isin = sum(k.(x).*w)
#     @test isapprox(I,Isin,rtol=1e-2)
# end

