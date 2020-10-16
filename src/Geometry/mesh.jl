"""
    abstract type AbstractMesh{N,T}
    
An abstract mesh structure in dimension `N` with primite data of type `T`. 

By default `T=Float64`.

An instance of an `AbstractMesh` is expected to implement the following methods:

#TODO: decide on the interface
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

function Base.eltype(M::AbstractMesh{N,T}) where {N,T}
    T
end

"""
    struct GenericMesh{N,T}

A simple data structure representing a generic mesh in an ambient space of dimension `N`, with data of type `T`. 
"""
struct GenericMesh{N,T} <: AbstractMesh{N,T}
    vtx::Vector{Point{N,T}}
    # element types
    etypes::Vector{Int32}
    # mapping from element type to indices of vtx in each element
    el2vtx::Vector{Matrix{Int}}
    # mapping from elementary entity to (etype,tags)
    ent2tags::Dict{ElementaryEntity,Dict{Int32,Vector{Int}}}
end

# TODO: write a method to `convert` a 3d mesh to a 2d mesh by projecting points
# into the firs two components

"""
    etypes(M::GenericMesh)

Return the element types contained in the mesh `M`.
"""
etypes(M::GenericMesh) = [type_tag_to_etype[i] for i in M.etypes]



"""
    const type_tag_to_etype

Dictionary mapping `gmsh` element types, given as `Int32`, to the internal
equivalent of those. 

Such a mapping is useful for generating function barriers in order to dispatch on
methods which work on a concrete subtype. 
"""
const type_tag_to_etype = Dict(
    15 => Point{3,Float64},
    1  => LagrangeLine{2},
    2  => LagrangeTriangle{3},
    4  => LagrangeTetrahedron{4}
)


