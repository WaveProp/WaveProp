"""
    struct TensorLagInterp{N,Td,T}

Generic Lagrange interpolation over an `N`-dimensional tensor grid. The
implementation uses a multidimensional generalization of the barycentric
formula.

The main constructor takes an `SVector{N,Vector{Td}}` containig the `N`
one-dimensional `nodes` and an `Array{N,T}` of the function `vals` at the tensor
product grid formed by the one-dimensional `nodes`.

# Examples:
```julia
nx = 10
ny = 12
x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
f   = (x) -> cos(x[1]*x[2])
vals = [f((x,y)) for x in x, y in y]
p   = TensorLagInterp(SVector(x,y),vals)
p((0.1,0.2)) ≈ f((0.1,0.2))
```
"""
struct TensorLagInterp{N,Td,T}
    vals::Array{T,N}
    nodes::SVector{N,Vector{Td}}
    weights::SVector{N,Vector{Td}}
end

nodes(p::TensorLagInterp)   = p.nodes
weights(p::TensorLagInterp) = p.weights
vals(p::TensorLagInterp)    = p.vals

return_type(::TensorLagInterp{_,__,T}) where {_,__,T} = T

ambient_dimension(::TensorLagInterp{N}) where {N} = N

function TensorLagInterp(nodes::SVector{N,Vector{Td}},vals::Array{T,N}) where {N,Td,T}
    weights = map(barycentric_lagrange_weights,nodes)
    TensorLagInterp{N,Td,T}(vals,nodes,weights)
end
TensorLagInterp(nodes::Vector,vals::Vector) = TensorLagInterp((nodes,),vals)
TensorLagInterp(nodes::NTuple,vals::Array)  = TensorLagInterp(SVector(nodes),vals)

function (p::TensorLagInterp{N,Td,T})(x::SVector) where {N,Td,T}
    num = zero(T)
    den = zero(Td)
    for I in CartesianIndices(p.vals)
        wi     = one(Td)
        x_m_xi = one(Td)
        for d in 1:N
            wi     *= p.weights[d][I[d]]
            x_m_xi *= p.nodes[d][I[d]] - x[d]
        end
        # FIXME: implicilty assumes that x is not one of the interpolation
        # nodes. Division by zero if that is the case.
        num += wi/x_m_xi*p.vals[I]
        den += wi/x_m_xi
    end
    num/den
end
(p::TensorLagInterp)(x::NTuple) = p(SVector(x))
(p::TensorLagInterp{1})(x::Number) = p((x,))

# multidimensional version of algorithm on page 504 of
# https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function barycentric_lagrange_weights(x::Vector)
    n = length(x) - 1
    w = similar(x)
    w[1] = 1.0
    for j in 1:n
        for k in 0:j-1
            w[k+1] = (x[k+1] - x[j+1]) * w[k+1]
        end
        w[j+1] = prod(0:j-1) do k
            x[j+1] - x[k+1]
        end
    end
    for j in 0:n
        w[j+1] = 1/w[j+1]
    end
    return w
end
