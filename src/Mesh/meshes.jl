"""
    abstract type AbstractMesh{N,T}
    
An abstract mesh structure in dimension `N` with primite data of type `T`. By
default `T=Float64`, and `N=3`.

The `AbstractMesh` interface expects the following methods to be implemented:

- `etypes(msh)` : return a `Vector{AbstractElement}` representing the element
  types composing the mesh.
- `elements(msh,etype)` : return an iterator over the mesh elements of type `etype::AbstractElement`
- `msh[E][i]` : return the *i-th* element of type `E`.
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

Data structure representing a generic mesh in an ambient space of dimension `N`, with data of type `T`. 
"""
Base.@kwdef struct GenericMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{SVector{N,T}} = Vector{SVector{N,T}}()
    # element types
    etypes::Vector{DataType} = Vector{DataType}()
    # for each element type (key), get the data required to reconstruct the
    # elements (value)
    elements::Dict{DataType,Any} = Dict{DataType,Any}()
    # mapping from elementary entity to (etype,tags)
    ent2tags::Dict{AbstractEntity,Dict{DataType,Vector{Int}}} = Dict{AbstractEntity,Dict{DataType,Vector{Int}}}()
end

# convert a mesh to 2d by ignoring third component. Note that this also requires
# converting various element types to their 2d counterpart.
function convert_to_2d(mesh::GenericMesh{3})
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
convert_to_2d(::Type{SVector{3,T}}) where {T} = SVector{2,T}

"""
    etypes(M::GenericMesh)

Return the element types contained in the mesh.
"""
etypes(mesh::GenericMesh) = mesh.etypes

# submesh structure


"""
    elements(mesh::AbstractMesh,E::Type)

Return an iterator for iterating over all elements of `mesh` of type `E`.
"""
elements(mesh::AbstractMesh, E::Type)   = ElementIterator(mesh,E)

# utilities for iterating over elements of a mesh

struct ElementIterator{E,M} <: AbstractVector{E}
    mesh::M
end  
ElementIterator{E}(mesh::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(mesh)
ElementIterator(mesh,E) = ElementIterator{E}(mesh)

Base.getindex(m::GenericMesh,E::Type{<:AbstractElement}) = ElementIterator(m,E)

Base.size(iter::ElementIterator) = (length(iter),)

mesh(iter::ElementIterator) = iter.mesh

function Base.getindex(iter::ElementIterator,I)
    @assert eltype(I) <: Integer
    map(i->iter[i],I)    
end    

function Base.length(iter::ElementIterator{<:LagrangeElement,<:GenericMesh})
    msh               = mesh(iter)    
    E                 = eltype(iter)    
    tags::Matrix{Int} = msh.elements[E]
    Np, Nel           = size(tags)
    return Nel
end    

function Base.getindex(iter::ElementIterator{<:LagrangeElement,<:GenericMesh},i::Int)
    E                   = eltype(iter)    
    mesh                = iter.mesh    
    tags::Matrix{Int}   = mesh.elements[E]
    node_tags           = view(tags,:,i)
    vtx                 = view(mesh.nodes,node_tags)
    el                  = E(vtx)    
    return el
end    

function Base.length(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E               = eltype(iter)
    mesh            = iter.mesh    
    tags::Vector{E} = mesh.elements[E]    
    return length(iter.mesh.elements[E])
end    

function Base.getindex(iter::ElementIterator{<:ParametricElement,<:GenericMesh},i::Int)
    E               = eltype(iter)
    mesh            = iter.mesh    
    els::Vector{E} = mesh.elements[E]    
    return els[i]
end    

"""
    eltindices(m::GenericMesh,Ω,E)

Compute the element indices `idxs` of the elements of type `E` composing `Ω`, so that `m[E][idxs]` gives all
the elements of type `E` meshing `Ω`.
"""
function eltindices(m::GenericMesh,Ω,E::DataType)
    idxs = Int[]    
    for ent in entities(Ω)
        tags = get(m.ent2tags[ent],E,Int[])
        append!(idxs,tags)
    end    
    return idxs
end    

"""
    eltindices(m::GenericMesh,Ω)

Return a `Dict` with keys being the `etypes` of `m`, and values being the
indices of the elements in `Ω` of type `E`. 
"""
function eltindices(m::GenericMesh,Ω)
    dict = Dict{DataType,Vector{Int}}()
    for E in etypes(m)
        tags = eltindices(m,Ω,E)    
        if !isempty(tags)
            dict[E] = tags
        end
    end
    return dict
end    

"""
    _qrule_for_mesh(m,p)

Given a mesh `m`, create a dictionary mapping each element type of `m` to an
appropriate quadrature rule of order `p` over that element type.

See also [`_qrule_for_reference_shape`](@ref)
"""
function _qrule_for_mesh(m,order)
    Dict(E=>Integration._qrule_for_element(E,order) for E in etypes(m))
end    