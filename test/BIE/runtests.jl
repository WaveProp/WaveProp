using SafeTestsets

@safetestset "Pde tests" begin include("pde_test.jl") end


# @testset "Operators" begin
#     pde  = Helmholtz(;dim=3,k=1)
#     G    = SingleLayerKernel(pde)
#     Ω, M = WaveProp.IO.gmsh_sphere()
#     qrule = Gauss{1}()
#     Q    = quadgen(M,qrule)
#     S    = IntegralOperator{ComplexF64}(G,Q,Q)
# end