using LinearAlgebra
using WaveProp
using WaveProp.Nystrom
using WaveProp.Geometry
using WaveProp.Integration
using WaveProp.Mesh
using WaveProp.ParametricSurfaces
using StaticArrays
using Random
using Plots
plotlyjs()
Random.seed!(1)

## general definitions
xout = SVector(3,3)
ops = (
    Laplace(;dim=2),
    Helmholtz(;dim=2,k=1.2),
    Elastostatic(;dim=2,μ=2,λ=3)
)
pde = ops[2] # chose your pde here
T    = Nystrom.default_density_eltype(pde)
c    = rand(T)

## place to store dofs and error
ndofs = []
ee = []

npatches = [2,4,8,16,32,64]
qorder   =  5
for n in npatches
    ## generate mesh
    clear_entities!()
    Ω   = ParametricSurfaces.Circle() |> Domain
    Γ   = boundary(Ω)
    M   = meshgen(Γ,(n,))
    mesh = NystromMesh(view(M,Γ),order=qorder)
    push!(ndofs,length(dofs(mesh)))
    ## create exact solution
    u    = (qnode) -> SingleLayerKernel(pde)(xout,qnode)*c
    dudn = (qnode) -> transpose(DoubleLayerKernel(pde)(xout,qnode))*c
    γ₀u   = γ₀(u,mesh)
    γ₁u   = γ₁(dudn,mesh)
    γ₀u_norm = norm(norm.(γ₀u,Inf),Inf)

    ## compute error
    S,D  = Nystrom.single_doublelayer_dim(pde,mesh)
    e1    = Nystrom.error_interior_green_identity(S,D,γ₀u,γ₁u)/γ₀u_norm
    push!(ee,norm(e1,Inf))
end

fig = plot(ndofs,ee,xscale=:log10,yscale=:log10,m=:o,label="error",lc=:black)
plot!(fig,xlabel="n",ylabel="error")
for p in 1:5
    c = ee[end]*ndofs[end]^p
    plot!(fig,ndofs,c ./ ndofs.^(p),xscale=:log10,yscale=:log10,label="h^$p",ls=:dash)
end
display(fig)
