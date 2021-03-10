"""
    struct Density{T,S} <: AbstractVector{T}

Discrete density with values `vals` on the qnodes `qnodes(surface::S)`.
"""
struct Density{V,S} <: AbstractVector{V}
    vals::Vector{V}
    surface::S
end
Base.size(σ::Density,args...)     = size(σ.vals,args...)
Base.getindex(σ::Density,args...) = getindex(σ.vals,args...)
Base.setindex(σ::Density,args...) = setindex(σ.vals,args...)

Density(etype::DataType,surf)  = Density(zeros(etype,length(surf)),surf)
Density(pde::AbstractPDE,surf) = Density(default_density_eltype(pde),surf)

Base.:-(σ::Density) = Density(-σ.vals,σ.surface)

function Base.:\(A::AbstractMatrix,σ::Density)
    @assert size(A,1) == size(A,2)
    y = A\σ.vals
    return Density(y,σ.surface)
end

function IterativeSolvers.gmres!(σ::Density,A,μ::Density,args...;kwargs...)
    gmres!(σ.vals,A,μ.vals,args...;kwargs...)
end
IterativeSolvers.gmres(A,μ::Density,args...;kwargs...) = gmres!(zero(μ),A,μ,args...;kwargs...)


function γ₀(f,X)
    vals = [f(x) for x in qnodes(X)]
    return Density(vals,X)
end

function γ₁(dfdn,X)
    vals = [dfdn(x,n) for (x,n) in zip(qnodes(X),qnormals(X))]
    return Density(vals,X)
end
