"""
    abstract type AbstractMesh{N,T}
    
An abstract mesh structure in dimension `N` with primite data of type `T`. By
default `T=Float64`, and `N=3`.

Subtypes of `AbstractMesh` are expected to implement the following methods:
- `etypes(msh)` : return a `Vector{AbstractElement}` representing the element
  types composing the mesh.
- `elements(msh,etype)`: return an iterator over the mesh elements of type `etype::AbstractElement`
- `nodes(msh)` : return the nodes of the mesh associated with its DOF.
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

# we define the geometric dimension of the mesh to be the largest of the geometric
# dimension of its entities. 
geometric_dimension(M::AbstractMesh) = maximum(x -> geometric_dimension(x), etypes(M))

Base.eltype(M::AbstractMesh{N,T}) where {N,T} = T

Base.length(mesh::AbstractMesh) = length(nodes(mesh))

"""
    struct GenericMesh{N,T} <: AbstractMesh{N,T}

A simple data structure representing a generic mesh in an ambient space of dimension `N`, with data of type `T`. 
"""
Base.@kwdef struct GenericMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{Point{N,T}} = Vector{Point{N,T}}()
    # element types
    etypes::Vector{DataType} = Vector{DataType}()
    # for each element type (key), get the data required to reconstruct the
    # elements (value)
    elements::Dict{DataType,Any} = Dict{DataType,Any}()
    # mapping from elementary entity to (etype,tags)
    ent2tags::Dict{ElementaryEntity,Dict{DataType,Vector{Int}}} = Dict{ElementaryEntity,Dict{DataType,Vector{Int}}}()
end

# convert a mesh to 2d by ignoring third component. Note that this also requires
# converting various element types to their 2d counterpart.
function GenericMesh{2}(mesh::GenericMesh{3})
    @assert all(x -> geometric_dimension(x) < 3, etypes(mesh)) 
    T = eltype(mesh)
    # create new dictionaries for elements and ent2tags with 2d elements as keys
    elements  = empty(mesh.elements)
    ent2tags  = empty(mesh.ent2tags)
    for (E, tags) in mesh.elements
        E2d = convert_to_2d(E)    
        elements[E2d] = tags
    end
    for (ent, dict) in mesh.ent2tags
        new_dict = empty(dict)    
        for (E, tags) in dict
            E2d = convert_to_2d(E)    
            new_dict[E2d] = tags
        end
        ent2tags[ent] = new_dict
    end    
    # construct new 2d mesh
    GenericMesh{2,T}(;
        nodes=[x[1:2] for x in nodes(mesh)],
        etypes=convert_to_2d.(etypes(mesh)),
        elements=elements,
        ent2tags=ent2tags
    )
end

nodes(m::GenericMesh)    = m.nodes
elements(m::GenericMesh) = m.elements
ent2tags(m::GenericMesh) = m.ent2tags

function convert_to_2d(::Type{LagrangeElement{R,N,SVector{3,T}}}) where {R,N,T} 
    LagrangeElement{R,N,SVector{2,T}}
end
convert_to_2d(::Type{Point{3,T}}) where {T} = Point{2,T}

"""
    etypes(M::GenericMesh)

Return the element types contained in the mesh.
"""
etypes(mesh::GenericMesh) = mesh.etypes

# submesh structure
struct SubMesh{N,T} <: AbstractMesh{N,T}
    mesh::GenericMesh{N,T}
    domain::Domain
end

Base.view(m::GenericMesh,Ω::Domain) = SubMesh(m,Ω)
Base.view(m::GenericMesh,ent::ElementaryEntity) = SubMesh(m,Domain(ent))

function etypes(submesh::SubMesh)
    Ω, M = submesh.domain, submesh.mesh    
    ee = DataType[]
    for ent in entities(Ω)
        dict = M.ent2tags[ent]
        append!(ee, keys(dict))
    end    
    return unique!(ee)
end   

struct ElementIterator{E,M}
    mesh::M
end    
ElementIterator{E}(mesh::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(mesh)
ElementIterator(mesh,E) = ElementIterator{E}(mesh)

"""
    elements(mesh::AbstractMesh,E::Type)

Return an iterator for iterating over all elements of `mesh` of type `E`.
"""
elements(mesh::AbstractMesh, E::Type)   = ElementIterator(mesh,E)

Base.eltype(::Type{ElementIterator{E,M}}) where {E,M} = E

# iterator for generic mesh (not associated with a domain)
function Base.length(iter::ElementIterator{<:LagrangeElement,<:GenericMesh})
    E                 = eltype(iter)    
    tags::Matrix{Int} = iter.mesh.elements[E]
    Np, Nel           = size(tags)
    return Nel
end    

function Base.length(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E        = eltype(iter)
    mesh     = iter.mesh    
    tags::Vector{E} = mesh.elements[E]    
    return length(iter.mesh.elements[E])
end    

Base.iterate(iter::ElementIterator{<:LagrangeElement,<:GenericMesh}) = iterate(iter,1)
function Base.iterate(iter::ElementIterator{<:LagrangeElement,<:GenericMesh},state)
    E                   = eltype(iter)    
    mesh                = iter.mesh    
    tags::Matrix{Int}   = mesh.elements[E]
    if state > length(iter)
        return nothing
    else    
        node_tags   = view(tags,:,state)
        vtx         = view(mesh.nodes,node_tags)
        el          = E(vtx)    
        return el, state + 1
    end
end    

Base.iterate(iter::ElementIterator{<:ParametricElement,<:GenericMesh}) = iterate(iter,1)
function Base.iterate(iter::ElementIterator{<:ParametricElement,<:GenericMesh}, state)
    E        = eltype(iter)
    mesh     = iter.mesh    
    tags::Vector{E} = mesh.elements[E]    
    iterate(tags,state)
end    

# iterator for submesh. Filter elements based on domain
function Base.length(iter::ElementIterator{<:LagrangeElement,<:SubMesh})
    submesh    = iter.mesh
    Ω, M        = submesh.domain, submesh.mesh
    E          = eltype(iter)    
    # loop over all entities and count the number of elements of type E
    Nel        = 0
    for ent in entities(Ω)
        dict   = M.ent2tags[ent] 
        v      = dict[E] # get element tags for type E
        Nel += length(v)    
    end
    return Nel
end    

function Base.length(iter::ElementIterator{<:ParametricElement,<:SubMesh})
    submesh    = iter.mesh
    Ω, M        = submesh.domain, submesh.mesh
    E          = eltype(iter)    
    # loop over all entities and count the number of elements of type E
    Nel        = 0
    for ent in entities(Ω)
        dict   = M.ent2tags[ent] 
        v      = dict[E] # get element tags for type E
        Nel += length(v)    
    end
    return Nel
end    

Base.iterate(iter::ElementIterator{<:AbstractElement,<:SubMesh}) = iterate(iter,(1, 1))
function Base.iterate(iter::ElementIterator{<:AbstractElement,<:SubMesh}, state)
    # extract some constant fields for convenience
    submesh   = iter.mesh    
    Ω, M      = submesh.domain, submesh.mesh
    E         = eltype(iter)    
    ents      = entities(Ω)
    # inner iteration over ent
    n, i          = state
    ent          = ents[n]
    inner_state  = _iterate(M, E, ent, i)
    if inner_state === nothing
        if n == length(ents)
            return nothing
        else
            return iterate(iter, (n + 1, 1))
        end
    else 
        el, i = inner_state
        return el, (n, i)
    end    
end  

# helper iterator which dispatches based on the type of E. A type stable method
# is obtained by *tagging* the mesh.elements[E] metadata with its expected type
function _iterate(mesh::GenericMesh,E::Type{<:LagrangeElement},ent::ElementaryEntity,i::Int=1)
    # if ent has not elements of type `E`, stop inner iteration
    haskey(mesh.ent2tags[ent],E) || (return nothing)
    # get tag for i-th elements in ent of type E
    el_tags = mesh.ent2tags[ent][E]
    if i > length(el_tags)
        return nothing    
    else
        el_tag      = el_tags[i]            
        tags::Matrix{Int} = mesh.elements[E]
        node_tags   = view(tags,:,el_tag)
        vtx         = view(mesh.nodes,node_tags)
        el          = E(vtx)
        return el,i+1
    end
end 

function _iterate(mesh::GenericMesh,E::Type{<:ParametricElement},ent::ElementaryEntity,i::Int=1)
    # if ent has not elements of type `E`, stop inner iteration
    haskey(mesh.ent2tags[ent],E) || (return nothing)
    # get tag for i-th elements in ent of type E
    el_tags = mesh.ent2tags[ent][E]
    if i > length(el_tags)
        return nothing    
    else
        el_tag      = el_tags[i]            
        tags::Vector{E} = mesh.elements[E]    
        el              = tags[el_tag]
        return el,i+1
    end
end 