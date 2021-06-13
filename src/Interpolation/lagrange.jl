# algorithm on page 504 of https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function barycentric_lagrange_weights(x)
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

# equation 4.2 of https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function barycentric_lagrange_matrix(xs,xt,w=barycentric_lagrange_weights(xs))
    ns = length(xs)
    nt = length(xt)
    d = similar(xt)
    for i in 1:nt
        d[i] = sum(zip(xs,w)) do (x,w)
            w /(xt[i]-x)
        end
    end
    A= [w[j]/(xt[i]-xs[j]) / d[i] for i in 1:nt, j in 1:ns]
    return A
end
