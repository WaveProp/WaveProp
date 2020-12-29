"""
    struct ElementIterator{E,M}

Return an iterator for iterating over all elements of type `E` on objects of
type `M`.
"""
struct ElementIterator{E,M} <: AbstractVector{E}
    mesh::M
end  

ElementIterator{E}(mesh::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(mesh)
ElementIterator(mesh,E) = ElementIterator{E}(mesh)

mesh(iter::ElementIterator) = iter.mesh

function Base.getindex(iter::ElementIterator,I)
    @assert eltype(I) <: Integer
    map(i->iter[i],I)    
end    

function Base.size(iter::ElementIterator{<:LagrangeElement,<:GenericMesh})
    msh               = mesh(iter)    
    E                 = eltype(iter)    
    tags::Matrix{Int} = msh.elements[E]
    Np, Nel           = size(tags)
    return (Nel,)
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

function Base.size(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E               = eltype(iter)
    mesh            = iter.mesh    
    tags::Vector{E} = mesh.elements[E]    
    return (length(iter.mesh.elements[E]),)
end    

function Base.getindex(iter::ElementIterator{<:ParametricElement,<:GenericMesh},i::Int)
    E               = eltype(iter)
    mesh            = iter.mesh    
    els::Vector{E} = mesh.elements[E]    
    return els[i]
end    

# iterator for submesh
function Base.size(iter::ElementIterator{<:AbstractElement,<:SubMesh})
    E          = eltype(iter)        
    submesh    = iter.mesh
    idxs       = dom2elt(submesh,E)
    return (length(idxs),)
end  

function Base.getindex(iter::ElementIterator{<:AbstractElement,<:SubMesh},i::Int)
    E      = eltype(iter)
    msh    = mesh(iter) # a SubMesh
    p_msh  = parent(msh) # parent mesh
    iglob  = dom2elt(msh,E)[i] # global index of element in parent mesh
    iter = ElementIterator(p_msh,E)
    return iter[iglob]
end    

ElementIterator(m::CartesianMesh) = ElementIterator(m,etype(m))

# element iterator for cartesian mesh
# FIXME: the ElementIterator for a CartesianMesh should inherit from
# AbstractMatrix and not AbstractVector since these are more naturally indexed
# as matrices. Also this should be made generic on the dimension
function Base.size(iter::ElementIterator{<:Any,<:CartesianMesh})
    E       = eltype(iter)    
    grids   = grid1d(iter.mesh)
    sz      = size(iter.mesh) 
    return sz .- 1
end    

# 1d case
function Base.getindex(iter::ElementIterator{<:Any,<:CartesianMesh{1}}, i::Int)
    E           = eltype(iter)
    xx          = xgrid(iter.mesh)
    low_corner  = (xx[i], )
    high_corner = (xx[i+1], )
    el          = E(low_corner,high_corner) # construct the element
    return el
end    

# 2d case
function Base.getindex(iter::ElementIterator{<:Any,<:CartesianMesh{2}}, i::Int,j::Int)
    E      = eltype(iter)    
    xx     = xgrid(iter.mesh)
    yy     = ygrid(iter.mesh)
    low_corner = (xx[i],yy[j])
    high_corner = (xx[i+1],yy[j+1])
    el  = E(low_corner,high_corner) # construct the element
    return el
end    

function Base.getindex(iter::ElementIterator{<:Any,<:CartesianMesh{2}}, I::Int)
    sz  = size(iter)
    idx = CartesianIndices(sz)[I]
    iter[idx[1],idx[2]]
end    