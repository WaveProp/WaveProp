struct SubMesh{N,T} <: AbstractMesh{N,T}
    mesh::GenericMesh{N,T}
    domain::Domain
    eltindices::Dict{DataType,Vector{Int}} # FIXME: this could be computed on the flight. Is that a problem?
    function SubMesh{N,T}(mesh::GenericMesh,Ω::Domain) where {N,T}
        idxs = eltindices(mesh,Ω)    
        return new{N,T}(mesh,Ω,idxs)
    end    
end
SubMesh(m::GenericMesh{N,T},args...;kwargs...) where {N,T} = SubMesh{N,T}(m,args...;kwargs...)

parent(m::SubMesh) = m.mesh

"""
    eltindices(m::SubMesh,[E])

A dictionary with keys being the element types of `m`, and values being the
element indices in the parent mesh. If a type `E` is given, return the values
associated with that key.
"""
eltindices(m::SubMesh) = m.eltindices
eltindices(m::SubMesh,E::Type{<:AbstractElement}) = m.eltindices[E]

Base.view(m::GenericMesh,Ω::Domain)           = SubMesh(m,Ω)
Base.view(m::GenericMesh,ent::AbstractEntity) = SubMesh(m,Domain(ent))

function etypes(submesh::SubMesh)
    Ω, M = submesh.domain, submesh.mesh    
    ee = DataType[]
    for ent in entities(Ω)
        dict = M.ent2tags[ent]
        append!(ee, keys(dict))
    end    
    return unique!(ee)
end   

# iterator for submesh
function Base.length(iter::ElementIterator{<:AbstractElement,<:SubMesh})
    E          = eltype(iter)        
    submesh    = iter.mesh
    idxs       = eltindices(submesh,E)
    return length(idxs)
end  

function Base.getindex(iter::ElementIterator{<:AbstractElement,<:SubMesh},i::Int)
    E      = eltype(iter)
    msh    = mesh(iter) # a SubMesh
    p_msh  = parent(msh) # parent mesh
    iglob  = eltindices(msh,E)[i] # global index of element in parent mesh
    return p_msh[E][iglob]
end    
