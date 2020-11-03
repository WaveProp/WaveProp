abstract type AbstractQuadrature{N,T} end

struct GenericQuadrature{N,T}
    nodes::Vector{Point{N,T}}
    weights::Vector{T}
    normals::Vector{Point{N,T}}
end

function GenericQuadrature{N,T}() where {N,T}
    GenericQuadrature{N,T}([],[],[])
end

qnodes(q::GenericQuadrature) = q.nodes
qweights(q::GenericQuadrature) = q.weights
qnormals(q::GenericQuadrature) = q.normals

function GenericQuadrature{N,T}(els::Tuple,qrule) where {N,T}
    q = GenericQuadrature{N,T}()
    for el in els
        x,w,ν = push_forward_quad_with_normal(el,qrule)
        append!(q.nodes,x)
        append!(q.weights,w)
        append!(q.normals,ν)
    end    
    return q
end    

struct SQuadrature{M,N,T}
    nodes::SVector{M,Point{N,T}}
    weights::SVector{M,T}
    normals::SVector{M,Point{N,T}}
end    
