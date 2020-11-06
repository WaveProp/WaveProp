using Test
using WaveProp.Nystrom

@testset "Pde basic tests" begin
    pde  = Helmholtz(;dim=3,k=1)
    @test Nystrom.ambient_dimension(pde) == 3
    pde  = Laplace(;dim=2)
    @test Nystrom.ambient_dimension(pde) == 2
end