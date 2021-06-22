"""
    struct Density{T,S} <: AbstractVector{T}

Discrete density with values `vals` on the quadrature nodes of `mesh::S`.
"""
struct Density{V,S<:NystromMesh} <: AbstractVector{V}
    vals::Vector{V}
    mesh::S
end

vals(σ::Density) = σ.vals
mesh(σ::Density) = σ.mesh

Base.size(σ::Density,args...)     = size(σ.vals,args...)
Base.getindex(σ::Density,args...) = getindex(σ.vals,args...)
Base.setindex(σ::Density,args...) = setindex(σ.vals,args...)

Density(etype::DataType,surf)  = Density(zeros(etype,length(surf)),surf)
Density(pde::AbstractPDE,surf) = Density(default_density_eltype(pde),surf)

Base.:-(σ::Density) = Density(-σ.vals,σ.mesh)

function Base.:\(A::AbstractMatrix,σ::Density)
    @assert size(A,1) == size(A,2)
    y = A\σ.vals
    return Density(y,σ.mesh)
end

function IterativeSolvers.gmres!(σ::Density,A,μ::Density,args...;kwargs...)
    gmres!(σ.vals,A,μ.vals,args...;kwargs...)
end
IterativeSolvers.gmres(A,μ::Density,args...;kwargs...) = gmres!(zero(μ),A,μ,args...;kwargs...)


function γ₀(f,X)
    vals = [f(dof) for dof in dofs(X)]
    return Density(vals,X)
end

function γ₁(dfdn,X)
    vals = [dfdn(dof) for dof in dofs(X)]
    return Density(vals,X)
end
