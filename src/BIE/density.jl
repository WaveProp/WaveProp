"""
    struct Density{T,S} <: AbstractVector{T}
    
Discrete density with values `vals` on the nodes `get_nodes(surface::S)`. 
"""
struct Density{T,S} <: AbstractVector{T}
    vals::Vector{T}
    surface::S
end    

Base.size(σ::Density) = size(σ.vals)
Base.getindex(σ::Density,args...) = getindex(σ.vals,args...)
Base.setindex(σ::Density,args...) = setindex(σ.vals,args...)

# zero density
Density{T}(X) where {T} = Density(zeros(T,length(X)),X)

function γ₀(f,X)
    vals = [f(x) for x in get_nodes(X)]
    isconcretetype(eltype(vals)) || warn("type inference failed: performance will be degraded")
    return Density(vals,X)
end

function γ₁(dfdn,X)
    vals = [dfdn(x,nx) for (x,nx) in zip(getnodes(X),getnormals(X))]
    isconcretetype(eltype(vals)) || warn("type inference failed: performance will be degraded")
    return Density(vals,X)
end