function elementary_matrix!(A,el::AbstractElement, u::PolynomialBasis,
                           v::PolynomialBasis, q::AbstractQuadratureRule;
                           f=(u,v)->(i,j,x̂,_,_)->u(x̂)[j]*v(x̂)[i])
    x̂,ŵ = q()
    x,w = q(el)
    m, n = length(v), length(u)
    for k in 1:length(w)
        for i in 1:m
            for j in 1:n
                A[i,j] += f(u, v)(i, j, x̂[k], el, x[k]) * w[k]
            end
        end
    end
    return A
end    

function elementary_matrix(el::AbstractElement, u::PolynomialBasis,
    v::PolynomialBasis, q::AbstractQuadratureRule;
    f=(u,v)->(i,j,x̂,_,_)->u(x̂)[j]*v(x̂)[i])
    x̂,ŵ = q()
    x,w = q(el)
    m, n = length(v), length(u)
    A = zeros(Float64,m,n)
    elementary_matrix!(A,el,u,v,q,f=f)
end

function elementary_matrix(el::AbstractElement, v::PolynomialBasis,
                           q::AbstractQuadratureRule;
                           f=(v)->(i,x̂,_,_)->v(x̂)[i])
    x̂,ŵ = q()
    x,w = q(el)
    m = length(v)
    b = zeros(Float64,m)
    for k in 1:length(w)
        for i in 1:m
            b[i] += f(v)(i, x̂[k], el, x[k]) * w[k]
        end
    end
    return b
end    

# NOTE: for testing purposes regarding speed differences between generic and a
# specific elementary matrix kernel
function (u::LagrangeBasis)(q::AbstractQuadratureRule)
    x̂,_ = q()
    mapreduce(x->u(x),hcat,x̂) |> transpose
end    

function mass_matrix(el::AbstractElement,u::PolynomialBasis,v::PolynomialBasis,q::AbstractQuadratureRule)
    x̂,ŵ = q()
    U = u(q)
    V = v(q)
    μ = map(u->measure(el,u),x̂)
    Σ = diagm(μ)*diagm(ŵ)
    return transpose(V)*Σ*U
end    

function stiffness_matrix!(A,el::AbstractElement,u::PolynomialBasis,v::PolynomialBasis,q::AbstractQuadratureRule)

end    
