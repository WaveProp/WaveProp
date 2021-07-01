using WaveProp
using WaveProp.Geometry
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using WaveProp.Nystrom
using WaveProp.Integration
using StaticArrays
using LinearAlgebra
using IterativeSolvers
using Random
using Plots

##
plotlyjs()
Random.seed!(1)
const k = 1.3
pde    = Maxwell(;dim=3,k)
qorder = 4

# exact solution
G    = SingleLayerKernel(pde)
xs   = SVector(0.1,0.2,0.3)
x0   = SVector(5,5,5)
c    = SVector(1+im,-2.,3.)
E    = (dof) -> G(dof,xs)*c
exa  = E(x0)

# geometry and parameters
order = 7
geo = ParametricSurfaces.Sphere(;radius=3)
Ω = Domain(geo)
Γ = boundary(Ω)
##

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
    er    = (Eₐ(x0) - exa)/norm(exa)
    @info norm(er,Inf)
end
