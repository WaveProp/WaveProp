function elementary_matrix(el,u::LagrangeBasis,v::LagrangeBasis,q::AbstractQuadratureRule)
    x̂,ŵ = q()
    x,w = q(el)
    U   = u(x̂)
    V   = v(x̂)
    return V*diagm(w)*transpose(U)
end    