using Test
using WaveProp.Nystrom
using StaticArrays

@testset "Kernels" begin
    pde  = Helmholtz(;dim=3,k=1)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == ComplexF64
    @test Nystrom.return_type(dG) == ComplexF64
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)

    pde  = Laplace(;dim=3)
    G    = SingleLayerKernel(pde)
    dG    = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == Float64
    @test Nystrom.return_type(dG) == Float64
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)

    pde  = Elastostatic(;μ=1,λ=2.0,dim=3)
    G    = SingleLayerKernel(pde)
    dG   = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == SMatrix{3,3,Float64,9}
    @test Nystrom.return_type(dG) == SMatrix{3,3,Float64,9}
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)

    pde  = Maxwell(;k=1.0)
    G    = SingleLayerKernel(pde)
    dG   = DoubleLayerKernel(pde)
    @test Nystrom.kernel_type(G) == Nystrom.SingleLayer()
    @test Nystrom.kernel_type(dG) == Nystrom.DoubleLayer()
    @test Nystrom.return_type(G) == SMatrix{3,3,ComplexF64,9}
    @test Nystrom.return_type(dG) == SMatrix{3,3,ComplexF64,9}
    @test Nystrom.combined_field_coefficients(G) == (0,-1)
    @test Nystrom.combined_field_coefficients(dG) == (1,0)
end

@testset "Maxwell identities" begin
    # essentially vector calculus identities verified by Maxwell's diadic greens
    # function
    op = Maxwell(2)
    T  = Nystrom.default_density_eltype(op)
    E  = rand(T)
    x  = rand(SVector{3})
    y  = rand(SVector{3})
    jac = rand(SMatrix{3,2})
    xdof = NystromDOF(x,1.,jac,0,0)
    ydof = NystromDOF(y,1.,jac,0,0)
    # [ n(y) × (∇ʸ × G(x,y) ) ]ᵗ E = [ ∇ʸ × G ]ᵗ (E × n(y))
    lhs = dG(xdof,ydof)*E
    rhs = Nystrom._curl_y_green_tensor_maxwell(x,y,op.k) * cross(E,normal(ydof))
    @test lhs ≈ rhs
    # [ ∇ʸ × G ]ᵗ = ∇ˣ × G
    lhs =  Nystrom._curl_y_green_tensor_maxwell(x,y,op.k) |> transpose
    rhs =  Nystrom._curl_x_green_tensor_maxwell(x,y,op.k)
    @test lhs == rhs
end
