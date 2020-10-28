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
geometric_dimension(M::AbstractMesh) = maximum(x->geometric_dimension(x),etypes(M))

Base.eltype(M::AbstractMesh{N,T}) where {N,T} = T

Base.length(mesh::AbstractMesh) = length(nodes(mesh))

"""
    struct GenericMesh{N,T} <: AbstractMesh{N,T}

A simple data structure representing a generic mesh in an ambient space of dimension `N`, with data of type `T`. 
"""
struct GenericMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{Point{N,T}}
    # element types
    etypes::Vector{DataType}
    # for each element type, the indices of nodes in each element
    el2nodes::Vector{Matrix{Int}}
    # mapping from elementary entity to (etype,tags)
    ent2tags::Dict{ElementaryEntity,Dict{Int32,Vector{Int}}}
    # quadrature info
    qnodes::Vector{Point{N,T}}
    qweights::Vector{T}
    qnormals::Vector{Point{N,T}}
    # for each element type, the indices of quadrature in each element
    el2qnodes::Vector{Matrix{Int}}
end

function GenericMesh(nodes, etypes, el2nodes, ent2tag)
    P = eltype(nodes)
    N,T = length(P), eltype(P)
    el2qnodes = Vector{Matrix{Int}}(undef,length(etypes))
    GenericMesh{N,T}(nodes,etypes,el2nodes,ent2tag,P[],T[],P[],el2qnodes)
end    

function GenericMesh{2}(mesh::GenericMesh{3})
    @assert all(x->geometric_dimension(x)<3,etypes(mesh)) 
    T = eltype(mesh)
    GenericMesh{2,T}(
        [x[1:2] for x in nodes(mesh)],
        convert_to_2d.(etypes(mesh)),
        mesh.el2nodes,
        mesh.ent2tags,
        [x[1:2] for x in qnodes(mesh)],
        qweights(mesh),
        [x[1:2] for x in qnormals(mesh)],
        mesh.el2qnodes
    )
end

convert_to_2d(::Type{LagrangeElement{R,N,3,T}}) where {R,N,T} = LagrangeElement{R,N,2,T}
convert_to_2d(::Type{Point{3,T}}) where {T} = Point{2,T}

nodes(m::GenericMesh)     = m.nodes
qnodes(m::GenericMesh)    = m.qnodes
qweights(m::GenericMesh)  = m.qweights
qnormals(m::GenericMesh)  = m.qnormals
el2qnodes(m::GenericMesh) = m.el2qnodes
el2qnodes(m::GenericMesh,i::Int) = m.el2qnodes[i]
function el2qnodes(m::GenericMesh,E::DataType) 
    i = findfirst(x->x==E,etypes(m))
    el2qnodes(m,i)
end

"""
    etypes(M::GenericMesh)

Return the element types contained in the mesh.
"""
etypes(M::GenericMesh) = M.etypes 

function elements(mesh::GenericMesh,E::Type{<:AbstractElement})
    return ElementIterator{E}(mesh)
end    

struct ElementIterator{E,M}
    mesh::M
end    
ElementIterator{E}(mesh::M) where {E,M<:AbstractMesh} = ElementIterator{E,M}(mesh)

Base.eltype(::Type{ElementIterator{E,M}}) where {E,M} = E

function Base.length(iter::ElementIterator{<:Any,<:GenericMesh})
    i       = _compute_etype_index(iter)
    tags    = iter.mesh.el2nodes[i]
    Np, Nel = size(tags)
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
    tags   = mesh.el2nodes[i]
    if state > length(iter)
        return nothing
    else    
        el_nodes = mesh.nodes[tags[:,state]] # get the coordinates of nodes in this element
        el  = E(el_nodes)                    # construct the element
        return el, state+1
    end
end    

"""
    _compute_quadrature!(msh::GenericMesh,E,qrule;need_normal=false)

For all elements of `msh` of type `E::Type{<:AbstractElement}`, use `qrule::AbstractQuadratureRule`
to compute an element quadrature and push that information into `msh`. Set
`need_normal=true` if the normal vector at the quadrature nodes should be computed.
"""
function _compute_quadrature!(mesh::GenericMesh,E,qrule;need_normal=false)
    @assert domain(qrule) == domain(E) "quadrature rule must be defined on domain of element"    
    N,T = ambient_dimension(mesh), eltype(mesh)
    i   = findfirst(x -> x==E, etypes(mesh))
    i === nothing && (return mesh)
    x̂,ŵ = qrule() # quadrature on reference element
    nq  = length(x̂) # number of qnodes per element
    el2qnodes = Int[]
    for el in elements(mesh,E)
        x,w = push_forward_quad(el,qrule)
        # compute indices of quadrature nodes in this element
        qidxs  = length(mesh.qnodes) .+ (1:nq) |> collect
        append!(el2qnodes,qidxs)
        append!(mesh.qnodes,x)
        append!(mesh.qweights,w)
        if need_normal==true
            n⃗ = map(u->normal(el,u),x̂) 
            append!(mesh.qnormals,n⃗)
        end
    end
    el2qnodes = reshape(el2qnodes, nq, :)
    mesh.el2qnodes[i] = el2qnodes
    return mesh   
end    

function _compute_quadrature!(mesh::GenericMesh,e2qrule::Dict;need_normal=false)
    for (E,qrule) in e2qrule
        _compute_quadrature!(mesh,E,qrule;need_normal)    
    end
    return mesh
end    

"""
    compute_quadrature!(mesh;order,dim,need_normal=false)

Compute a quadrature of a desired `order` for all elements of dimension `dim`.
Set `need_normal=true` if the normal at the quadrature nodes is to be computed.
"""
function compute_quadrature!(mesh::GenericMesh;order,dim,need_normal=false)
    dict = Dict()
    @assert allunique(etypes(mesh))
    for E in etypes(mesh)
        geometric_dimension(E) == dim || continue
        ref = domain(E)
        qrule = _qrule_for_reference_element(ref,order)
        push!(dict,E=>qrule)
    end
    _compute_quadrature!(mesh,dict;need_normal)
    return mesh
end 

"""
    _qrule_for_reference_element(ref,order)

Given a `ref`erence element and a desired quadrature `order`, return 
an appropiate quadrature rule.
"""
function _qrule_for_reference_element(ref,order)
    if ref isa ReferenceLine
        n = ((order + 1) ÷  2) + 1
        qrule = GaussLegendre{n}()
    elseif ref isa ReferenceSquare
        n  = (order + 1)/2 |> ceil
        qx = GaussLegendre{n}()
        qy = qx
        qrule = TensorProduct(qx,qy)
    elseif ref isa ReferenceTriangle
        if order <= 1
            return Gauss(ref,n=1) 
        elseif order <=2
            return Gauss(ref,n=3)     
        else
            @notimplemented    
        end
    elseif ref isa ReferenceTetrahedron
        if order <= 1
            return Gauss(ref;n=1) 
        elseif order <=2
            return Gauss(ref;n=4)
        else
            @notimplemented    
        end
    end    
end    






