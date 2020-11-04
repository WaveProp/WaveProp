using Test
using QuadGK
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.SingularIntegration
using LinearAlgebra

@testset "QuadGK" begin
    function test_quadgk(f,sing_handler,a,b,xs)
        val,er  = quadgk(f,a,xs,b)
        val2,er = quadgk(sing_handler,f,a,xs,b,xs=xs)
        @test val ≈ val2
    end
    f = (x) -> x==0.1 ? 0.0 : log(abs(x-0.1))
    kress = Kress{10}()
    test_quadgk(f,kress,-1,1,0.1)
    imt = IMT{1,1}()
    test_quadgk(f,imt,-1,1,0.1)
end

@testset "Duffy" begin
    @testset "Singularity at vertex" begin
        q1d   = GaussLegendre{5}()
        qstd  = TensorProductQuadrature(q1d,q1d)
        duffy = Duffy{2}()
        qsin = SingularQuadratureRule(qstd,duffy)
        x,w  = qsin()
        # regular kernel
        k(x) = cos(x[1])
        I    = quadgk(k,ReferenceTriangle)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-5)
        # singular kernel at right vertex
        s    = Point(1,0) 
        k(x) = 1/norm(x-s)
        I    = quadgk(k,ReferenceTriangle)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-5)
        # singular kernel at wrong vertex (should not be accurate, see test below)
        s    = Point(0,1) 
        k(x) = 1/norm(x-s)
        I    = quadgk(k,ReferenceTriangle)
        Isin = sum(k.(x).*w)
        @test !isapprox(I,Isin,rtol=1e-5)
    end

    @testset "Singularity inside" begin
        q1d   = GaussLegendre{5}()
        qstd  = TensorProductQuadrature(q1d,q1d)
        duffy = Duffy{2}()
        qsin = SingularQuadratureRule(qstd,duffy)
        s    = Point(0.37,0.25) 
        # triangle case
        x,w  = qsin(ReferenceTriangle(),s)
        # regular kernel
        k(x) = cos(x[1])
        I    = quadgk(k,ReferenceTriangle)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-5)
        # singular kernel at interior point
        k(x) = 1/norm(x-s)
        I    = quadgk(k,ReferenceTriangle)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-2)
        # square case
        x,w  = qsin(ReferenceSquare(),s)
        # regular kernel
        k(x) = cos(x[1])
        I    = quadgk(k,ReferenceSquare)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-5)
        # singular kernel at interior point
        k(x) = 1/norm(x-s)
        I    = quadgk(k,ReferenceSquare)
        Isin = sum(k.(x).*w)
        @test isapprox(I,Isin,rtol=1e-2)
    end    
end

@testset "2d Kress" begin
    q1d   = GaussLegendre{10}()
    qstd  = TensorProductQuadrature(q1d,q1d)
    s1d   = Kress(order=2)
    sing_handler = TensorProductQuadratureHandler(s1d,s1d)
    qsin = SingularQuadratureRule(qstd,sing_handler)
    x,w  = qsin()    
    # regular kernel
    k(x) = cos(x[1])
    I    = quadgk(k,ReferenceSquare)
    Isin = sum(k.(x).*w)
    @test isapprox(I,Isin,rtol=1e-2)
    # singular kernel at origin
    k(x) = 1/norm(x)
    I    = quadgk(k,ReferenceSquare)
    Isin = sum(k.(x).*w)
    @test isapprox(I,Isin,rtol=1e-2)
    # singular kernel elsewhere
    s   = Point(0.47,0.38)
    x,w  = qsin(ReferenceSquare(),s)
    k(x) = 1/norm(x-s)
    I    = quadgk(k,ReferenceSquare)
    Isin = sum(k.(x).*w)
    @test isapprox(I,Isin,rtol=1e-2)
end

