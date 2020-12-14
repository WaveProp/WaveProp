using Test
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using LinearAlgebra

@testset "Lagrange elements" begin
    @testset "LagrangeLine" begin
        ## line in 3d
        vtx = ((0.,0.,0.),(1.,1.,1.))
        l   = LagrangeLine(vtx)
        @test domain(l) == ReferenceLine()
        @test ambient_dimension(l) == 3
        @test geometric_dimension(l) == 1
        # makes no sense to ask for normal here. Make sure error is thrown.
        @test_throws AssertionError normal(l,SVector(0.5)) == SVector(-1.,1.)/√2
        ## line in 2d
        a = (0.,0.)
        b = (1.,1.)
        l   = LagrangeLine((a,b))
        @test domain(l) == ReferenceLine()
        @test ambient_dimension(l) == 2
        @test geometric_dimension(l) == 1
        @test normal(l,SVector(0.5)) == SVector(1.,-1.)/√2
    end
    @testset "LagrangeTriangle" begin
        # triangle in 2d
        vtx = SVector(SVector(0.,0.),SVector(0.,1.),SVector(-1.,0))
        t   = LagrangeTriangle(vtx)
        @test ambient_dimension(t) == 2
        @test geometric_dimension(t) == 2
        # triangle in 3d
        vtx = SVector(SVector(0.,0.,0.),SVector(0.,1.,0.),SVector(-1.,0,0.))
        t   = LagrangeTriangle(vtx)
        @test ambient_dimension(t) == 3
        @test geometric_dimension(t) == 2
        @test normal(t,SVector(0.1,0.1)) == SVector(0,0,1.)
    end
    @testset "Tetrahedron" begin
        # TODO: add tests
    end
end

@testset "Line quadrature" begin
    el     = LagrangeLine((1.,1.),(5.,4.))
    qrule  = GaussLegendre{1}()
    x,w    = qrule(el)
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

@testset "Curved lines" begin
    qstd = GaussLegendre(20)    
    el   = LagrangeLine((0,0),(1,1),(1/2,1/4))
    # simple test on smooth integrand
    f     = x->cos(x[1])
    I     = integrate(f,el)
    for shand in [IMT(), Kress()]
        q     = SingularQuadratureRule(qstd,shand)
        Ia    = integrate(f,q,el)
        @debug Ia - I
        @test isapprox(Ia,I,rtol=1e-6)
    end
    # non-smooth integrand
    f     = x->log(abs(x[1])*cos(x[1]))
    I     = integrate(f,el)
    for shand in [IMT(), Kress()]
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
    for shand in [IMT(), Kress()]
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