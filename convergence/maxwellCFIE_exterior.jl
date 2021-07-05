using WaveProp
using WaveProp.Geometry
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using WaveProp.Nystrom
using WaveProp.Integration
using StaticArrays
using LinearAlgebra
using IterativeSolvers
using Plots
plotlyjs()

# Convergence test using Carlos'
# operators and definitions

##
const k = 1.3
pde = MaxwellCFIE(k)

# geometry and parameters
geo = ParametricSurfaces.Sphere(;radius=3)
Ω = Domain(geo)
Γ = boundary(Ω)

# evaluation mesh
eval_geo = ParametricSurfaces.Sphere(;radius=5)
eval_Ω = Domain(eval_geo)
eval_Γ = boundary(eval_Ω)
eval_M = meshgen(eval_Γ,(2,2))
eval_nystrom_mesh = NystromMesh(view(eval_M,eval_Γ);order=4)
eval_mesh = qcoords(eval_nystrom_mesh) |> collect

# exact solution
G    = (x,y) -> Nystrom.maxwell_green_tensor(x, y, k)
curlG = (x,y) -> Nystrom.maxwell_curl_green_tensor(x, y, k)
xs   = SVector(0.1,0.2,0.3)
c    = SVector(1+im,-2.,3.)
E    = (dof) -> G(dof,xs)*c
curlE = (dof) -> curlG(dof,xs)*c
exa  = E.(eval_mesh);

## Direct formulation
for n in [6]
    M     = meshgen(Γ,(n,n))
    mesh  = NystromMesh(view(M,Γ);order=5)
    γ₀E   = ncross(trace(E,mesh))
    γ₁E   = ncross(trace(curlE,mesh))
    S,D   = Nystrom.single_doublelayer_dim(pde,mesh;n_src=26)
    L     = I/2 - D
    rhs   = S*γ₁E
    γ₀Eₐ   = L\rhs
    γ₁Eₐ   = γ₁E
    Spot  = Nystrom.maxwellCFIE_SingleLayerPotencial(pde, mesh)
    Dpot  = Nystrom.maxwellCFIE_DoubleLayerPotencial(pde, mesh)
    Eₐ    = (x) -> Spot(γ₁Eₐ, x) + Dpot(γ₀Eₐ, x)
    er    = (Eₐ.(eval_mesh) - exa)/norm(exa,Inf)
    @info norm(er,Inf)
    #ϕ_gmres   =  zero(rhs)
    #gmres!(ϕ_gmres,L,rhs;verbose=true,maxiter=600,restart=600,abstol=1e-6)
end

## Exterior Green identity
u    = (qnode) -> SingleLayerKernel(pde)(qnode,xs)*c  # n × G
dudn = (qnode) -> DoubleLayerKernel(pde)(qnode,xs)*c  # n × ∇ × G
for n in [1,2,4,6]
    M     = meshgen(Γ,(n,n))
    mesh  = NystromMesh(view(M,Γ);order=5)
    γ₀u   = trace(u,mesh)
    γ₁u   = trace(dudn,mesh)
    γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)
    S, D = Nystrom.single_doublelayer_dim(pde,mesh;n_src=26)
    e1 = (S*γ₁u + D*γ₀u - 0.5*γ₀u) / γ₀u_norm  # exterior Green identity
    @info norm(e1,Inf)
end

## Indirect formulation
n = 4
M     = meshgen(Γ,(n,n))
mesh  = NystromMesh(view(M,Γ);order=5)
γ₀E   = ncross(trace(E,mesh))
S,D   = Nystrom.single_doublelayer_dim(pde,mesh;n_src=26)
_,J = WaveProp.Nystrom.diagonal_ncross_jac_matrix(mesh)
L     = Matrix{SMatrix{2,2,ComplexF16,4}}(transpose(J)*(I/2 + D)*J)
rhs   = transpose(J)*γ₀E
ϕ_coeff =  L\rhs
ϕ     = J*ϕ_coeff
Dpot  = Nystrom.maxwellCFIE_DoubleLayerPotencial(pde, mesh)
Eₐ    = (x) -> Dpot(ϕ,x)
er    = (Eₐ.(eval_mesh) - exa)/norm(exa,Inf)
@info norm(er,Inf)

#ϕ_gmres   =  zero(rhs)
#gmres!(ϕ_gmres,L,rhs;verbose=true,maxiter=600,restart=600,abstol=1e-6)


