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
import WaveProp.Nystrom: coords
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
exa  = E.(eval_mesh)

## Direct formulation
for n in [1,2,4,6]
    M     = meshgen(Γ,(n,n))
    mesh  = NystromMesh(view(M,Γ);order=5)
    γ₀E   = trace(E,mesh)
    S,D   = Nystrom.single_doublelayer_dim(pde,mesh)
    γ₁H   = -im*k*ncross(γ₀E)
    L     = -I/2 + D
    rhs   = S*γ₁H
    γ₀H   =  zero(rhs)
    # gmres!(γ₀H,L,rhs;verbose=true,maxiter=600,restart=600)
    γ₀H   =  L\rhs
    γ₁Eₐ  = im*k*ncross(γ₀H)
    γ₀Eₐ  = ncross(γ₁H)/(im*k)
    Spot  = SingleLayerPotential(pde,mesh)
    Dpot  = DoubleLayerPotential(pde,mesh)
    Eₐ    = (x) -> Dpot[γ₀Eₐ](x) - Spot[γ₁Eₐ](x)
    er    = (Eₐ.(eval_mesh) - exa)/norm(exa,Inf)
    @info norm(er,Inf)
end
