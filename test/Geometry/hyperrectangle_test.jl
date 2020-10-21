using Test
using WaveProp.Geometry

@testset "HyperRectangle tests" begin
    low_corner  = (0.0,0.0)
    high_corner = (1.0,2.0)
    mid         = (low_corner .+ high_corner) ./ 2 |> Point
    rec = HyperRectangle(low_corner,high_corner)
    @test mid == center(rec)
    @test (mid ∈ rec)  == true
    @test !in(high_corner .+ (1,1),rec)
    @test diameter(rec) == sqrt(1^2 + 2^2)
    @test radius(rec) == sqrt(1^2 + 2^2)/2
    # bbox
    pts = NTuple{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,(x,y))
        end
    end
    @test bounding_box(pts) == HyperRectangle((-1.,-1),(1,1.))
end
