"""
    struct SubMesh{N,T} <: AbstractMesh{N,T}

Create a view of a `parent` mesh over a given `domain`.

A submesh implements the interface for `AbstractMesh`; therefore you can iterate
over elements of the submesh just like you would with a mesh.
"""
struct SubMesh{N,T} <: AbstractMesh{N,T}
    parent::GenericMesh{N,T}
    domain::Domain
    dom2elt::Dict{DataType,Vector{Int}}
    function SubMesh{N,T}(mesh::GenericMesh,Ω::Domain) where {N,T}
        idxs = dom2elt(mesh,Ω)
        return new{N,T}(mesh,Ω,idxs)
    end
end
SubMesh(m::GenericMesh{N,T},args...;kwargs...) where {N,T} = SubMesh{N,T}(m,args...;kwargs...)

Base.view(m::GenericMesh,Ω::Domain)           = SubMesh(m,Ω)
Base.view(m::GenericMesh,ent::AbstractEntity) = SubMesh(m,Domain(ent))

# TODO: write tests for this
Base.view(m::SubMesh,Ω::Domain)           = view(mesh(m),intersect(Ω,domain(m)))
Base.view(m::SubMesh,ent::AbstractEntity) = view(m,Domain(ent))

parent(m::SubMesh) = m.parent
domain(m::SubMesh) = m.domain

# ElementIterator for submesh
function Base.size(iter::ElementIterator{<:AbstractElement,<:SubMesh})
    E          = eltype(iter)
    submesh    = mesh(iter)
    idxs       = dom2elt(submesh,E)
    return (length(idxs),)
end

function Base.getindex(iter::ElementIterator{<:AbstractElement,<:SubMesh},i::Int)
    E      = eltype(iter)
    submsh = mesh(iter) # a SubMesh
    p_msh  = parent(submsh) # parent mesh
    iglob  = dom2elt(submsh,E)[i] # global index of element in parent mesh
    iter   = ElementIterator(p_msh,E) # iterator over parent mesh
    return iter[iglob]
end

function Base.iterate(iter::ElementIterator{<:AbstractElement,<:SubMesh}, state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end

"""
    dom2elt(m::SubMesh,[E])

A dictionary with keys being the element types of `m`, and values being the
element indices in the parent mesh. If a type `E` is given, return the values
associated with that key.
"""
dom2elt(m::SubMesh) = m.dom2elt
dom2elt(m::SubMesh,E::Type{<:AbstractElement}) = m.dom2elt[E]

function Base.keys(submesh::SubMesh)
    Ω, M = submesh.domain, submesh.parent
    ee = DataType[]
    for ent in entities(Ω)
        dict = M.ent2tags[ent]
        append!(ee, keys(dict))
    end
    return unique!(ee)
end
