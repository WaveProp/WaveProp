function elementary_matrix(el::AbstractElement, u::PolynomialBasis,
                           v::PolynomialBasis, q::AbstractQuadratureRule;
                           f=(u,v)->(i,j,x̂,_,_)->u(x̂)[j]*v(x̂)[i])
    x̂,ŵ = q()
    x,w = q(el)
    m, n = length(v), length(u)
    A = zeros(Float64,m,n)
    for k in 1:length(w)
        for i in 1:m
            for j in 1:n
                A[i,j] += f(u, v)(i, j, x̂[k], el, x[k]) * w[k]
            end
        end
    end
    return A
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