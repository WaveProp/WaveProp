# Single Layer kernel for Maxwell is the dyadic Greens function
function (SL::SingleLayerKernel{T,S})(target,source)::T where {T,S<:MaxwellCFIE}
    k = parameters(pde(SL))
    x  = coords(target)
    y  = coords(source)
    rvec = x - y
    r = norm(rvec)
    r == 0 && return zero(T)
    # helmholtz greens function
    g   = exp(im*k*r)/(4π*r)
    gp  = im*k*g - g/r
    gpp = im*k*gp - gp/r + g/r^2
    RRT = rvec*transpose(rvec) # rvec ⊗ rvecᵗ
    G   = g*I + 1/k^2*(gp/r*I + (gpp/r^2 - gp/r^3)*RRT)
    # TODO: when multiplying by a density, it is faster to exploit the outer
    # product format isntead of actually assemblying the matrices.
    nx = cross_product_matrix(normal(target))
    return nx*G
end 

# Double Layer Kernel
function (DL::DoubleLayerKernel{T,S})(target,source)::T where {T,S<:MaxwellCFIE}
    k = parameters(pde(DL))
    x  = coords(target)
    y  = coords(source)
    rvec = x - y
    r      = norm(rvec)
    r == 0 && return zero(T)
    g      = exp(im*k*r)/(4π*r)
    gp     = im*k*g - g/r
    rcross = cross_product_matrix(rvec)
    nx = cross_product_matrix(normal(target))
    return nx * rcross * gp/r
end

function maxwell_green_tensor(target, source, k)
    x  = coords(target)
    y  = coords(source)
    rvec = x - y
    r = norm(rvec)
    r == 0 && return zero(T)
    # helmholtz greens function
    g   = exp(im*k*r)/(4π*r)
    gp  = im*k*g - g/r
    gpp = im*k*gp - gp/r + g/r^2
    RRT = rvec*transpose(rvec) # rvec ⊗ rvecᵗ
    G   = g*I + 1/k^2*(gp/r*I + (gpp/r^2 - gp/r^3)*RRT)
    return G
end

function maxwell_curl_green_tensor(target, source, k)
    x  = coords(target)
    y  = coords(source)
    rvec = x - y
    r      = norm(rvec)
    r == 0 && return zero(T)
    g      = exp(im*k*r)/(4π*r)
    gp     = im*k*g - g/r
    rcross = cross_product_matrix(rvec)
    return rcross * gp/r
end

function maxwellCFIE_SingleLayerPotencial(pde, mesh)
    k = parameters(pde)
    Γ = mesh
    function out(σ, x)
        iter = zip(dofs(Γ),σ)
        return sum(iter) do (source,σ)
            w = weight(source)
            maxwell_green_tensor(x, source, k)*σ*w
        end
    end
    return out
end
function maxwellCFIE_DoubleLayerPotencial(pde, mesh)
    k = parameters(pde)
    Γ = mesh
    function out(σ, x)
        iter = zip(dofs(Γ),σ)
        return sum(iter) do (source,σ)
            w = weight(source)
            maxwell_curl_green_tensor(x, source, k)*σ*w
        end
    end
    return out
end