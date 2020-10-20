"""
    quadgen(el::AbstractElement,qrule::AbstractQuadratureRule) -> (x,w)

Generate a quadrature for `el` using `qrule`. Return the quadrature nodes `x` and weights `w`.

The quadrature is computed by generating a quadrature on the reference
element of `el`, followed by the push-forward map. This requires:

    - el(x)
    - jacobian(el,x)
"""
function quadgen(el::AbstractElement,qrule)
    # generate a quadrature on the reference element    
    msg = "the domain of `qrule` must coincide with the domain `el`"
    @assert Geometry.domain(el) == domain(qrule) msg
    x̂,ŵ = qrule()
    # modify the quadrature using the push-forward map
    x   = map(el,x̂)
    μ   = map(x̂) do x̂
        jac = jacobian(el,x̂)
        g   = transpose(jac)*jac |> det
        sqrt(g)
    end 
    w   = ŵ .* μ 
    return x,w
end    
# FIXME: the function above is somewhat inneficient when the ambient and
# geometric dimensions of the element are the same. In that case `μ` simplifies
# to the usual `|det(jac)|`. This should be easy to fix by checking e.g. whether
# `jac` is a square matrix. Since these are static arrays there should be no
# runtime overhead compared to the hand-written version

function quadgen(M::GenericMesh,qrule;dim=ambient_dimension(M),need_normal=false)
    N,T = ambient_dimension(M), eltype(M)
    Q = GenericQuadrature{N,T}()
    for (i,etype) in enumerate(etypes(M))
        geometric_dimension(etype) == dim || continue 
        tags     = M.el2vtx[i]  # Np × Nel matrix
        Np, Nel  = size(tags)   # num of pts per element, num. of elements
        for n in 1:Nel
            el_vtx = M.vtx[tags[:,n]] # get the coordinates of nodes in this element
            el  = etype(el_vtx)       # construct the element
            @assert Geometry.domain(el) == domain(qrule)
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