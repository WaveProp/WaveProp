using Test
using WaveProp.BIE

@testset "Pde basic tests" begin
    pde  = Helmholtz(;dim=3,k=1)
    @test BIE.ambient_dimension(pde) == 3
    pde  = Laplace(;dim=2)
    @test BIE.ambient_dimension(pde) == 2
end