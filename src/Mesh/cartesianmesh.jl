"""
    struct UniformCartesianMesh{N,T} <: AbstractMesh{N,T}

An `N`-dimensional cartesian grid given as the tensor-product of `N`
one-dimensional `LinRange{T}` grids.

Iterating over a `UniformCartesianMesh` generates the elements which compose the mesh;
i.e. the `HyperRectangle` cells.
"""
struct UniformCartesianMesh{N,T} <: AbstractMesh{N,T}
    grids::NTuple{N,LinRange{T}}
end

grids(g::UniformCartesianMesh)     = g.grids
grids(g::UniformCartesianMesh,dim) = g.grids[dim]

Base.keys(m::UniformCartesianMesh{N,T}) where {N,T} = (HyperRectangle{N,T},)

"""
    UniformCartesianMesh(;domain::HyperRectangle,sz::NTuple)

Construct a uniform `UniformCartesianMesh` with `sz[d]` elements along dimension `d`.
"""
function UniformCartesianMesh(domain::HyperRectangle{N,T},sz::NTuple{N}) where {N,T}
    lc = low_corner(domain)
    uc = high_corner(domain)
    grids1d = ntuple(N) do n
        # to have sz elements, need sz+1 points in LinRange.
        npts = sz[n] + 1
        LinRange{float(T)}(lc[n], uc[n], npts) # use float(T) in case T<:Integer
    end
    UniformCartesianMesh(grids1d)
end
UniformCartesianMesh(domain::HyperRectangle{N,T},sz::Int) where {N,T} = UniformCartesianMesh(domain,ntuple(i->sz,N))
UniformCartesianMesh(;domain,sz) = UniformCartesianMesh(domain,sz)

ambient_dimension(g::UniformCartesianMesh{N}) where {N} = N

xgrid(g::UniformCartesianMesh) = g.grids[1]

ygrid(g::UniformCartesianMesh) = g.grids[2]

zgrid(g::UniformCartesianMesh) = g.grids[3]

# implement ElementIterator interface to UniformCartesianMesh

function Base.size(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh})
    g = mesh(iter)
    length.(g.grids) .- 1
end

Base.length(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh})
    CartesianIndices(size(iter))
end

function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    lc = ntuple(N) do dim
        i = I[dim]
        m.grids[dim][i]
    end
    hc = ntuple(N) do dim
        i = I[dim] + 1
        m.grids[dim][i]
    end
    HyperRectangle(lc, hc)
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh},I...)
    iter[CartesianIndex(I)]
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    iter[I[i]]
end

function Base.iterate(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh},state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end

# since UniformCartesianMesh has only one elment type, calling ElementIterator without
# specifying the has clear sense
function ElementIterator(m::UniformCartesianMesh)
    E = keys(m) |> first
    ElementIterator(m,E)
end

# NodeIterator over UniformCartesianMesh

function Base.IteratorSize(iter::Type{NodeIterator{UniformCartesianMesh{N,T}}}) where {N,T}
    Base.HasShape{N}()
end

function Base.size(iter::NodeIterator{<:UniformCartesianMesh})
    g = mesh(iter)
    length.(g.grids)
end

Base.length(iter::NodeIterator{<:UniformCartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::NodeIterator{<:UniformCartesianMesh})
    CartesianIndices(size(iter))
end

function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    svector(N) do dim
        i = I[dim]
        m.grids[dim][i]
    end
end
function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh},I...)
    iter[CartesianIndex(I)]
end
function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    iter[I[i]]
end

function Base.iterate(iter::NodeIterator{<:UniformCartesianMesh},state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end
