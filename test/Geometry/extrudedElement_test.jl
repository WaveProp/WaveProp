using Test
using WaveProp.Geometry
using WaveProp.Integration
using LinearAlgebra

@testset "Extruded line" begin    
    l = line(Point(1,1),Point(2.,2))
    d⃗ = normal(l,0.5)
    translate(l,d⃗)
    Ω = extrude(l,d⃗)
    @test geometric_dimension(l) == 1
    @test geometric_dimension(Ω) == 2
    ∂Ω    = boundary(Ω)
    qrule = GaussLegendre{10}()
    x,w = quadgen(∂Ω[1],qrule)
    @test sum(w) ≈  √2
    x,w = quadgen(∂Ω[2],qrule)
    @test sum(w) ≈  norm(d⃗)
    x,w = quadgen(∂Ω[3],qrule)
    @test sum(w) ≈  √2
    x,w = quadgen(∂Ω[4],qrule)
    @test sum(w) ≈ norm(d⃗)
end