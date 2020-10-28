using Test
using WaveProp.Geometry
using WaveProp.Utils

@testset "IMT" begin
    l   = ReferenceLine()
    n   = 10
    x̂   = svector(i->Point(i/n),n)
    cov = IMT()
    x   = map(cov,x̂)
    @test map(x-> x ∈ l,x) |> all
    μ   = map(x->jacobian(cov,x),x̂)
end

@testset "Duffy" begin
    s   = ReferenceSquare()
    t   = ReferenceTriangle()
    n,m   = 10,10
    x̂   = [Point((i-1)/(n-1),(j-1)/(m-1)) for i in 1:n, j in 1:m]
    x̂   = SVector{n*m}(x̂)
    @test map(x-> x ∈ s,x̂) |> all
    cov = Duffy{2}()
    x   = map(cov,x̂)
    @test map(x-> x ∈ t,x) |> all
end