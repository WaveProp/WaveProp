"""
    abstract type AbstractMesh{N,T}

An abstract mesh structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

Concrete subtypes of `AbstractMesh` should implement [`ElementIterator`](@ref)
for accessing the mesh elements, and/or [`NodeIterator`](@ref) for accesing the
mesh nodes.

# See also: [`GenericMesh`](@ref), [`UniformCartesianMesh`](@ref)
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

# we define the geometric dimension of the parent to be the largest of the geometric
# dimension of its entities.
geometric_dimension(M::AbstractMesh) = maximum(x -> geometric_dimension(x), entities(M))

primitive_type(M::AbstractMesh{N,T}) where {N,T} = T

function Base.show(io::IO,msh::AbstractMesh)
    print(io,"$(typeof(msh)) containg:")
    for E in keys(msh)
        iter = ElementIterator(msh,E)
        print(io,"\n\t $(length(iter)) elements of type ",E,)
    end
    return io
end

"""
    struct ElementIterator{E,M}

Return an iterator for all elements of type `E` on a mesh of type `M`.

Besides the methods listed in the [iterator
iterface](https://docs.julialang.org/en/v1/manual/interfaces/) of `Julia`, some
functions also require the `getindex(iter,i::Int)` method for accessing the
`i`-th element directly.
"""
struct ElementIterator{E,M}
    mesh::M
end

mesh(iter::ElementIterator) = iter.mesh

Base.eltype(::SType{ElementIterator{E}}) where {E} = E

ElementIterator{E}(parent::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(parent)

ElementIterator(parent,E) = ElementIterator{E}(parent)

# indexing a mesh with an element type return an interator over that element
# (like a dict)
Base.getindex(m::AbstractMesh,E::DataType) = ElementIterator(m, E)

function Base.length(iter::ElementIterator)
    prod(size(iter))
end

"""
    struct NodeIterator{M}

Iterator for all the nodes in a mesh of type `M`.
"""
struct NodeIterator{M}
    mesh::M
end

mesh(iter::NodeIterator) = iter.mesh

# nodes of `AbstractMesh{N,T}` should be points of type `SVector{N,T}`
Base.eltype(::SType{NodeIterator{<:AbstractMesh{N,T}}}) where {N,T} = SVector{N,T}
