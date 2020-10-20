"""
    struct Density{T,S} <: AbstractVector{T}
    
Discrete density with values `vals` on the nodes `get_nodes(surface::S)`. 
"""
struct Density{V,S} <: AbstractVector{V}
    vals::Vector{V}
    surface::S
end
Base.size(σ::Density,args...)     = size(σ.vals,args...)
Base.getindex(σ::Density,args...) = getindex(σ.vals,args...)
Base.setindex(σ::Density,args...) = setindex(σ.vals,args...)

Density(etype::DataType,surf) = Density(zeros(etype,length(surf)),surf)

function IterativeSolvers.gmres!(σ::Density,A,μ::Density,args...;kwargs...)
    gmres!(σ.vals,A,μ.vals,args...;kwargs...)
    return σ
end

function γ₀(f,X)
    vals = [f(x) for x in nodes(X)]
    return Density(vals,X)
end

function γ₁(dfdn,X)
    vals = [dfdn(x,n) for (x,n) in zip(nodes(X),normals(X))]
    return Density(vals,X)
end
