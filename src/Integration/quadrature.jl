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
    x̂,ŵ = qrule()
    for el in els
        x,w,n = _push_forward_quad_with_normal(el,x̂,ŵ)
        append!(q.nodes,x)
        append!(q.weights,w)
        append!(q.normals,n)
    end    
    return q
end    

function _push_forward_quad(el,x̂,ŵ)
    x   = map(x->el(x),x̂)
    w   = map(zip(x̂,ŵ)) do (x̂,ŵ)
        jac = jacobian(el,x̂)
        g   = transpose(jac)*jac |> det
        sqrt(g)*prod(ŵ)
    end 
    return x,w
end    

function _push_forward_quad_with_normal(el,x̂,ŵ)
    x,w = _push_forward_quad(el,x̂,ŵ)
    n   = map(x->normal(el,x),x̂)
    return x,w,n
end    
