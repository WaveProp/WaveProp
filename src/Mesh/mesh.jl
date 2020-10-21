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

Base.eltype(M::AbstractMesh{N,T}) where {N,T} = T

Base.length(mesh::AbstractMesh) = length(nodes(mesh))

"""
    struct GenericMesh{N,T} <: AbstractMesh{N,T}

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

"""
    etypes(M::GenericMesh)

Return the element types contained in the mesh.
"""
etypes(M::GenericMesh) = [type_tag_to_etype[i] for i in M.etypes]

function elements(mesh::GenericMesh,E::Type{<:AbstractElement})
    return ElementIterator{E}(mesh)
end    

struct ElementIterator{E,M}
    mesh::M
end    
ElementIterator{E}(mesh::M) where {E,M} = ElementIterator{E,M}(mesh)

# Base.eltype(iter::ElementIterator{E}) where {E} = E
Base.eltype(::Type{ElementIterator{E,M}}) where {E,M} = E

function Base.length(iter::ElementIterator{<:Any,<:GenericMesh})
    i       = _compute_etype_index(iter)
    tags    = iter.mesh.el2vtx[i]
    Np, Nel = size(tags)            # num of pts per element, num. of elements
    return Nel
end    

function _compute_etype_index(iter::ElementIterator{<:Any,<:Any})
    E = eltype(iter)
    i = findfirst(x -> x==E, etypes(iter.mesh))
    return i
end

function Base.iterate(iter::ElementIterator{<:Any,<:GenericMesh},state=1)
    E      = eltype(iter)    
    mesh   = iter.mesh    
    i      = _compute_etype_index(iter)
    tags   = mesh.el2vtx[i]
    if state > length(iter)
        return nothing
    else    
        el_vtx = mesh.vtx[tags[:,state]] # get the coordinates of nodes in this element
        el  = E(el_vtx)                  # construct the element
        return el, state+1
    end
end    

"""
    quadgen(M::GenericMesh,qrule;dim=ambient_dimension(M),need_normal=false)

Generate a quadrature using `qrule` for all elements of `M` having dimension
`dim`. The flag `need_normal` controls whether the normal at each quadrature
node is computed and stored.
"""
function quadgen(mesh::GenericMesh,qrule;dim=ambient_dimension(mesh),need_normal=false)
    N,T = ambient_dimension(mesh), eltype(mesh)
    Q = GenericQuadrature{N,T}()
    for (i,E) in enumerate(etypes(mesh))
        geometric_dimension(E) == dim || continue 
        for el in elements(mesh,E)
            @assert domain(el) == domain(qrule)
            x̂,ŵ = qrule()
            x,w = push_forward_map(el,x̂,ŵ)
            append!(Q.nodes,x)
            append!(Q.weights,w)
            if need_normal==true
                n⃗ = map(u->normal(el,u),x̂) 
                append!(Q.normals,n⃗)
            end
        end    
    end     
    return Q   
end    

struct NystromMesh{N,T} <: AbstractMesh{N,T}
    elements::Vector{AbstractElement}
    quadrature::GenericQuadrature{N,T}
    el2quad::Vector{Vector{Int}} # map from element to idx of quadrature nodes for that element
end    

function NystromMesh(mesh::GenericMesh,qrule;dim=ambient_dimension(mesh)-1)
    N,T = ambient_dimension(mesh), eltype(mesh)
    Q = GenericQuadrature{N,T}()
    els = Vector{AbstractElement}()
    el2quad = Vector{Vector{Int}}()
    for (i,etype) in enumerate(etypes(mesh))
        geometric_dimension(etype) == dim || continue 
        tags     = mesh.el2vtx[i]  # Np × Nel matrix
        Np, Nel  = size(tags)   # num of pts per element, num. of elements
        for n in 1:Nel
            el_vtx = mesh.vtx[tags[:,n]] # get the coordinates of nodes in this element
            el  = etype(el_vtx)       # construct the element
            push!(els,el)
            @assert domain(el) == domain(qrule)
            x̂,ŵ = qrule()
            x,w = push_forward_map(el,x̂,ŵ)
            nquad = length(x)
            idxs  = length(Q.nodes) .+ (collect(1:nquad))
            push!(el2quad,idxs)
            append!(Q.nodes,x)
            append!(Q.weights,w)
            n⃗ = map(u->normal(el,u),x̂) 
            append!(Q.normals,n⃗)
        end    
    end 
    return NystromMesh{N,T}(els,Q,el2quad)
end    

nodes(m::NystromMesh)    = m.quadrature.nodes
normals(m::NystromMesh)  = m.quadrature.normals
weights(m::NystromMesh)  = m.quadrature.weights
elements(m::NystromMesh) = m.elements