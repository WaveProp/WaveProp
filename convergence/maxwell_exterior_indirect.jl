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

##
const k = 1.3
pde = Maxwell(;dim=3,k)

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
G    = SingleLayerKernel(pde)
xs   = SVector(0.1,0.2,0.3)
c    = SVector(1+im,-2.,3.)
E    = (dof) -> G(dof,xs)*c
exa  = E.(eval_mesh);

## Indirect formulation (tangential density)
n = 4
M     = meshgen(Γ,(n,n))
mesh  = NystromMesh(view(M,Γ);order=5)
γ₀E   = trace(E,mesh)
S,D   = Nystrom.single_doublelayer_dim(pde,mesh)
N,J = WaveProp.Nystrom.diagonal_ncross_jac_matrix(mesh)
L     = Matrix{SMatrix{2,2,ComplexF16,4}}(transpose(J)*N*(I/2 + D)*J)
rhs   = transpose(J)*(N*γ₀E)
ϕ_coeff = L\rhs
ϕ     = J*ϕ_coeff
Dpot  = DoubleLayerPotential(pde,mesh)
Eₐ    = (x) -> Dpot[ϕ](x)
er    = (Eₐ.(eval_mesh) - exa)/norm(exa,Inf)
@info norm(er,Inf)

#ϕ_gmres   =  zero(rhs)
#gmres!(ϕ_gmres,L,rhs;verbose=true,maxiter=600,restart=600,abstol=1e-6)

## Indirect formulation
n = 4
M     = meshgen(Γ,(n,n))
mesh  = NystromMesh(view(M,Γ);order=5)
γ₀E   = trace(E,mesh)
S,D   = Nystrom.single_doublelayer_dim(pde,mesh)
N,_ = WaveProp.Nystrom.diagonal_ncross_jac_matrix(mesh)
L     = (I/2 + D)
rhs   = γ₀E
ϕ     = L\rhs
Dpot  = DoubleLayerPotential(pde,mesh)
Eₐ    = (x) -> Dpot[ϕ](x)
er    = (Eₐ.(eval_mesh) - exa)/norm(exa,Inf)
@info norm(er,Inf)

