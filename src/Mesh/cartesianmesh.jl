"""
    struct CartesianMesh{N,T} <: AbstractMesh{N,T}

An `N`-dimensional cartesian grid given as the tensor-product of `N` one-dimensional
`LinRange{T}` objects. The grid spacing is therefore constant per dimension.
"""
struct CartesianMesh{N,T} <: AbstractMesh{N,T}
    grid1d::NTuple{N,LinRange{T}}
end
etype(m::CartesianMesh{N,T}) where {N,T}  = HyperRectangle{N,T}
etypes(m::CartesianMesh) = (etype(m),)

# allow for arguments to be passed as e.g. (x,y) instead of tuple and promote
# type if needed
CartesianMesh(args...) = CartesianMesh(promote(args...))

function CartesianMesh(d::ReferenceLine;h::Number=floatmax(),n::Int=1)
    m     = cld(1,h) |> Int
    n     = max(m,n)
    xgrid = LinRange(0,1,n+1)
    CartesianMesh(xgrid)
end

function CartesianMesh(d::ReferenceSquare, h, n)
    nx     = cld(1,h[1]) + 1 |> Int
    ny     = cld(1,h[2]) + 1 |> Int
    nx     = max(n[1],nx)
    ny     = max(n[2],ny)
    xgrid = LinRange(0,1,nx)
    ygrid = LinRange(0,1,ny)
    CartesianMesh(xgrid,ygrid)
end
CartesianMesh(d::ReferenceSquare; h::Number, n::Number) = CartesianMesh(d,(h,h),(n,n))

grid1d(g::CartesianMesh)     = g.grid1d
grid1d(g::CartesianMesh,dim) = g.grid1d[dim]

dimension(g::CartesianMesh{N}) where {N} = N

xgrid(g::CartesianMesh) = g.grid1d[1]
ygrid(g::CartesianMesh) = g.grid1d[2]
zgrid(g::CartesianMesh) = g.grid1d[3]

meshsize(g::CartesianMesh)      = step.(grid1d(g))
meshsize(g::CartesianMesh,dim)  = step(grid1d(g,dim))

Base.size(g::CartesianMesh) = length.(g.grid1d)
Base.length(g::CartesianMesh) = prod(size(g))

function Base.getindex(g::CartesianMesh,I)
    N = dimension(g)
    @assert N == length(I)
    ntuple(N) do dim
        i = I[dim]
        g.grid1d[dim][i]
    end
end

function Base.getindex(g::CartesianMesh,I...)
    N = dimension(g)
    @assert N == length(I)
    ntuple(N) do dim
        i = I[dim]
        g.grid1d[dim][i]
    end
end

Base.CartesianIndices(g::CartesianMesh) = CartesianIndices(size(g))

# iterate over all nodes
# function Base.iterate(g::CartesianMesh)
#     i = first(CartesianIndices(g))
#     return g[i],i
# end

# function Base.iterate(g::CartesianMesh,state)
#     idxs = CartesianIndices(g)
#     next = iterate(idxs,state)
#     if next === nothing
#         return nothing
#     else
#         i,state = next
#         return g[i],state
#     end
# end

# Base.IteratorSize(::CartesianMesh{N}) where {N} = Base.HasShape{N}()
