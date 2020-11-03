using Test
using QuadGK
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.SingularIntegration
using WaveProp.Interpolation

@testset "Singular quadrature rules" begin
    qstd = GaussLegendre{10}()
    cov  = Kress()
    qsin = SingularQuadratureRule(qstd,cov)
    s = 0.1
    # create a singular function 
    f = (x) -> log(abs(x[1]-s))*cos(x[1])
    I,E = quadgk(f,0,s,1)
    Istd = integrate(f,qstd)
    Isin = integrate(f,qsin,s)
    e_std = I - Istd |> abs
    e_sin = I - Isin |> abs    
    @test !isapprox(I,Istd,rtol=1e-5)
    @test isapprox(I,Isin,rtol=1e-5)
    k = (x) -> log(abs(x[1]-s))
    f = (x) -> cos(x[1])
    I,E = quadgk(x->k(x)*f(x),0,s,1)
    x,w = qsin(k,s)
    Isin = sum(f.(x).*w)
    @test isapprox(Isin,I,rtol=1e-5)
end

@testset "Singular weights reference space" begin
    qstd = GaussLegendre{10}()
    cov  = Kress()
    qsin = SingularQuadratureRule(GaussLegendre{10}(),cov)
    s = 0.1
    k = (x) -> log(abs(x[1]-s))
    f = (x) -> cos(x[1])
    xi = qstd()[1]
    w  = singular_weights(qsin,xi,k,s)
    I,E = quadgk(x->k(x)*f(x),0,s,1)
    fi = f.(xi)
    Isin = sum(fi .* w )
    @test isapprox(Isin,I,rtol=1e-5)
end

