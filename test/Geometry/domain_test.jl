using Test
using WaveProp.Geometry

# Two points
p1 = ElementaryEntity(0, 1, ElementaryEntity[])
p2 = ElementaryEntity(0, 2, ElementaryEntity[])
points = [p1, p2]
# Three lines
l1 = ElementaryEntity(1, 1, ElementaryEntity[p1, p2])
l2 = ElementaryEntity(1, 2, ElementaryEntity[p1, p2])
l3 = ElementaryEntity(1, 3, ElementaryEntity[p1, p2])
lines = [l1, l2, l3]
# Two surfaces
s1 = ElementaryEntity(2, 1, ElementaryEntity[l1, l2])
s2 = ElementaryEntity(2, 2, ElementaryEntity[l2, l3])
surfaces = [s1, s2]

@testset "ElementaryEntity" begin
    @test tag(p1) == (0,1)
    @test geometric_dimension(p1) == 0
    @test boundary(l1) == points
    for e in vcat(points, lines, surfaces)
        for b in boundary(e)
            @test geometric_dimension(e)-1 == geometric_dimension(b)
        end
    end
end

# Domains
Ω1 = Domain(s1)
Ω2 = Domain(s2)
Ω = Domain(surfaces)

@testset "Domain" begin
    @test entities(Ω) == surfaces
    @test length(Ω) == 2
    @test !isempty(Ω)
    @test Ω2 == setdiff(Ω, Ω1)
    @test Ω == union(Ω1, Ω2)
    @test s1 in Ω
    @test Ω[1] == s1
    @test Ω[end] == s2
    for ω in Ω
        for Γ in boundary(ω)
            @test geometric_dimension(ω)-1 == geometric_dimension(Γ)
        end
    end
    @test Ω1 == intersect(Ω1, Ω)
    @test Domain(l2) == intersect(Ω1, Ω2)
    @test issubset(Ω1, Ω)
    @test Ω2 == remove(Ω1, Ω)
    @test union(Domain.(lines)...) == skeleton(Ω)
    @test Domain(l2) == internal_boundary(Ω)
    @test union(Domain(l1), Domain(l3)) == external_boundary(Ω)
    @test tags(Ω, 1) == [(1,1), (1,2), (1,3)]
    @test tags(Ω) == [(2,1), (2,2)]
end