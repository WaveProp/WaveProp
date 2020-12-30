struct SubMesh{N,T} <: AbstractMesh{N,T}
    mesh::GenericMesh{N,T}
    domain::Domain
    dom2elt::OrderedDict{DataType,Vector{Int}} # FIXME: this could be computed on the flight. Is that a problem?
    function SubMesh{N,T}(mesh::GenericMesh,Ω::Domain) where {N,T}
        idxs = dom2elt(mesh,Ω)    
        return new{N,T}(mesh,Ω,idxs)
    end    
end
SubMesh(m::GenericMesh{N,T},args...;kwargs...) where {N,T} = SubMesh{N,T}(m,args...;kwargs...)

Base.view(m::GenericMesh,Ω::Domain)           = SubMesh(m,Ω)
Base.view(m::GenericMesh,ent::AbstractEntity) = SubMesh(m,Domain(ent))

Base.getindex(m::GenericMesh,Ω::Domain) = view(m,Ω)
Base.getindex(m::GenericMesh,ent::AbstractEntity) = view(m,ent)

Base.getindex(m::SubMesh,Ω::Domain)           = view(parent(m),Ω)
Base.getindex(m::SubMesh,ent::AbstractEntity) = view(parent(m),ent)

parent(m::SubMesh) = m.mesh
domain(m::SubMesh) = m.domain

"""
    dom2elt(m::SubMesh,[E])

A dictionary with keys being the element types of `m`, and values being the
element indices in the parent mesh. If a type `E` is given, return the values
associated with that key.
"""
dom2elt(m::SubMesh) = m.dom2elt
dom2elt(m::SubMesh,E::Type{<:AbstractElement}) = m.dom2elt[E]

function etypes(submesh::SubMesh)
    Ω, M = submesh.domain, submesh.mesh    
    ee = DataType[]
    for ent in entities(Ω)
        dict = M.ent2tags[ent]
        append!(ee, keys(dict))
    end    
    return unique!(ee)
end   

