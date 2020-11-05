struct SingularQuadratureRule{Q,S}
    qrule::Q
    singularity_handler::S
end 

qrule(q::SingularQuadratureRule) = q.qrule
singularity_handler(q::SingularQuadratureRule) = q.singularity_handler

domain(q::SingularQuadratureRule) = q |> qrule |> domain
range(q::SingularQuadratureRule)  = q |> singularity_handler |> range

function (qs::SingularQuadratureRule)()
    qstd = qs |> qrule
    cov  = qs |> singularity_handler    
    x,w  = qstd(cov)
    return x,w
end

"""
    (qrule::SingularQuadratureRule)(k,s)

Return the nodes and weights to integrate `∫k(x)f(x)dx ≈ sum(f.(x) .* w)`.
"""
function (qrule::SingularQuadratureRule)(k,s)
    x,w = qrule(s)
    w   = map(zip(x,w)) do (x,w)
        k(x)*w
    end
    return x,w
end 

"""
    (qrule::SingularQuadratureRule)(s)

Return a quadrature rule to integrate a function over `domain(qrule)`. The
function can be (weakly) singular at the location `s`.
"""
function (qrule::SingularQuadratureRule)(s)
    x,w   = qrule()    
    # split the domain
    l1    = line(s,0)
    l2    = line(s,1)
    # apply the quadrature
    x1,w1 = qrule(l1)
    x2,w2 = qrule(l2)
    return vcat(x1,x2), vcat(w1,w2)
end 

function (qrule::SingularQuadratureRule{<:Any,Duffy})(::ReferenceTriangle,s)
    x,w   = qrule()    
    # split the domain
    t1    = triangle((0,0),s,(1,0))
    t2    = triangle((0,0),s,(0,1))
    t3    = triangle((1,0),s,(0,1))
    # apply the quadrature
    x1,w1 = qrule(t1)
    x2,w2 = qrule(t2)
    x3,w3 = qrule(t3)
    return vcat(x1,x2,x3), vcat(w1,w2,w3)
end 

function (qrule::SingularQuadratureRule{<:Any,Duffy})(::ReferenceSquare,s)
    x,w   = qrule()    
    # split the domain
    t1    = triangle((0,0),s,(1,0))
    t2    = triangle((0,0),s,(0,1))
    t3    = triangle((1,0),s,(1,1))
    t4    = triangle((0,1),s,(1,1))
    # apply the quadrature
    x1,w1 = qrule(t1)
    x2,w2 = qrule(t2)
    x3,w3 = qrule(t3)
    x4,w4 = qrule(t4)
    return vcat(x1,x2,x3,x4), vcat(w1,w2,w3,w4)
end 



"""
    singular_weights(q::SingularQuadratureRule,xi,k,s)

Return the weights to integrate `∫k(x)f(x)dx ≈ sum(f.(xi) .* w)` where a
singularity can be efficiently handled at location `s`.
"""
function singular_weights(qsin::SingularQuadratureRule,xi,k,s)
    x,w = qsin(k,s)
    wlag = barycentric_lagrange_weights(xi)
    L   = barycentric_lagrange_matrix(xi,x,wlag)
    return transpose(L)*w
end 

function singular_weights(qsin::SingularQuadratureRule,xi,K,s,el)
    k   = (v) -> begin
        jac = jacobian(el,v)
        μ   = transpose(jac)*jac |> det |> sqrt
        K(el(v))*μ
    end    
    singular_weights(qsin,xi,k,s)
end 

function integrate(f,q::SingularQuadratureRule)
    x,w = q()
    integrate(f,x,w)
end

function integrate(f,q::SingularQuadratureRule,s)
    x,w = q(s)
    integrate(f,x,w)
end    

function (qrule::SingularQuadratureRule{<:Any,<:TensorProductSingularityHandler})(::ReferenceSquare,s)
    sx,sy = s[1],s[2]
    # split the domain
    s1    = rectangle(s,(0,sy),(0,0),(sx,0))
    s2    = rectangle(s,(sx,0),(1,0),(1,sy))
    s3    = rectangle(s,(1,sy),(1,1),(sx,1))
    s4    = rectangle(s,(sx,1),(0,1),(0,sy))
    # apply the quadrature
    x1,w1 = qrule(s1)
    x2,w2 = qrule(s2)
    x3,w3 = qrule(s3)
    x4,w4 = qrule(s4)
    return vcat(x1,x2,x3,x4), vcat(w1,w2,w3,w4)
end 


