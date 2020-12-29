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

function CartesianMesh(d::ReferenceLine;gridsize::Number)
    n     = cld(1,gridsize) + 1 |> Int
    xgrid = LinRange(0,1,n)
    CartesianMesh(xgrid)
end    

function CartesianMesh(d::ReferenceSquare;gridsize)
    nx     = cld(1,gridsize[1]) + 1 |> Int
    ny     = cld(1,gridsize[2]) + 1 |> Int
    xgrid = LinRange(0,1,nx)
    ygrid = LinRange(0,1,ny)
    CartesianMesh(xgrid,ygrid)
end    

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
function Base.iterate(g::CartesianMesh)
    i = first(CartesianIndices(g))
    return g[i],i
end    

function Base.iterate(g::CartesianMesh,state)
    idxs = CartesianIndices(g)        
    next = iterate(idxs,state)            
    if next === nothing
        return nothing
    else    
        i,state = next
        return g[i],state
    end
end    

Base.IteratorSize(::CartesianMesh{N}) where {N} = Base.HasShape{N}()
