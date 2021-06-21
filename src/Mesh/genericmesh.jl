"""
    struct GenericMesh{N,T} <: AbstractMesh{N,T}

Data structure representing a generic mesh in an ambient space of dimension `N`,
with data of type `T`.
"""
Base.@kwdef struct GenericMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{SVector{N,T}} = Vector{SVector{N,T}}()
    # for each element type (key), store the data required to reconstruct the
    # elements (value).
    elements::Dict{DataType,Any} = Dict{DataType,Any}()
    # mapping from entity to (etype,tags)
    ent2tags::Dict{AbstractEntity,Dict{DataType,Vector{Int}}} = Dict{AbstractEntity,Dict{DataType,Vector{Int}}}()
end

nodes(m::GenericMesh)    = m.nodes
elements(m::GenericMesh) = m.elements
ent2tags(m::GenericMesh) = m.ent2tags

Base.keys(m::GenericMesh) = keys(elements(m))
entities(m::GenericMesh) = keys(ent2tags(m)) |> collect
domain(m::GenericMesh) = entities(m) |> Domain

# implement the interface for ElementIterator of lagrange elements on a generic mesh
function Base.size(iter::ElementIterator{<:LagrangeElement,<:GenericMesh})
    msh               = mesh(iter)
    E                 = eltype(iter)
    tags::Matrix{Int} = msh.elements[E]
    _, Nel           = size(tags)
    return (Nel,)
end

function Base.getindex(iter::ElementIterator{<:LagrangeElement,<:GenericMesh},i::Int)
    E                   = eltype(iter)
    M                   = mesh(iter)
    tags::Matrix{Int}   = M.elements[E]
    node_tags           = view(tags,:,i)
    vtx                 = view(M.nodes,node_tags)
    el                  = E(vtx)
    return el
end

function Base.iterate(iter::ElementIterator{<:LagrangeElement,<:GenericMesh}, state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end

"""
    dom2elt(m::GenericMesh,Ω,E)

Compute the element indices `idxs` of the elements of type `E` composing `Ω`, so that `m[E][idxs]` gives all
the elements of type `E` meshing `Ω`.
"""
function dom2elt(m::GenericMesh, Ω, E::DataType)
    idxs = Int[]
    for ent in entities(Ω)
        tags = get(m.ent2tags[ent], E, Int[])
        append!(idxs, tags)
    end
    return idxs
end

"""
    dom2elt(m::GenericMesh,Ω)

Return a `Dict` with keys being the element types of `m`, and values being the
indices of the elements in `Ω` of that type.
"""
function dom2elt(m::GenericMesh, Ω)
    dict = Dict{DataType,Vector{Int}}()
    for E in keys(m)
        tags = dom2elt(m, Ω, E)
        if !isempty(tags)
            dict[E] = tags
        end
    end
    return dict
end
