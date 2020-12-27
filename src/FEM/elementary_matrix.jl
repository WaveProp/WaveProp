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

function elementary_matrix_new!(A,el::AbstractElement, u::PolynomialBasis,
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
    mapreduce(x->u(x),hcat,x̂)
end    

function mass_matrix(el::AbstractElement,u::PolynomialBasis,v::PolynomialBasis,q::AbstractQuadratureRule)
    x̂,ŵ = q()
    U = u(q)
    V = v(q)
    μ = map(u->measure(el,u),x̂)
    Σ = diagm(μ)*diagm(ŵ)
    return V*Σ*transpose(U)
end    

function mass_matrix_unrolled!(A,el::AbstractElement,u::PolynomialBasis,v::PolynomialBasis,q::AbstractQuadratureRule)
    x̂,ŵ = q()
    U = u(q)
    V = v(q)
    U *= diagm(ŵ) # put precomputed weights here. Should be computed only once.
    for k in 1:length(ŵ)
        J = @inbounds jacobian(el,x̂[k])
        μ = abs(det(J))
        for j in 1:size(U,2)
            @inbounds uj = U[j,k]*μ    
            for i in 1:size(V,2)
                @inbounds A[i,j] += uj*V[i,k]
            end    
        end
    end    
    return A
end    

function (gu::gradLagrangeBasis)(q::AbstractQuadratureRule)
    x̂,_ = q()
    mapreduce(x->gu(x),hcat,x̂) 
end    

function stiffness_matrix_unrolled!(A,el::AbstractElement,u::PolynomialBasis,v::PolynomialBasis,q::AbstractQuadratureRule)
    x̂,ŵ = q()
    ∇U = grad(u)(q)
    ∇V = grad(v)(q)
    ∇U *= diagm(ŵ) # put the reference weights in one of the precomputed matrices
    for k in 1:length(ŵ)
        J  = jacobian(el,x̂[k])
        iJ = inv(J)
        μ  = abs(det(J))
        Σ  = iJ*transpose(iJ)*μ
        for i in 1:size(∇V,2)
            for j in 1:size(∇U,2)    
                A[i,j] += transpose(∇V[i,k]) * Σ * ∇U[j,k]
            end    
        end
    end    
    return A
end    
