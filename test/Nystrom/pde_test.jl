using Test
using WaveProp.Nystrom
using WaveProp.Geometry

@testset "Pde basic tests" begin
    pde  = Helmholtz(;dim=3,k=1)
    @test ambient_dimension(pde) == 3
    pde  = Laplace(;dim=2)
    @test ambient_dimension(pde) == 2
end