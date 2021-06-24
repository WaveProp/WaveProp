# Single Layer kernel for Maxwell is the dyadic Greens function
function (SL::SingleLayerKernel{T,S})(target,source)::T  where {T,S<:Maxwell}
    k  = parameters(pde(SL))
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
    return G
end

# Double Layer Kernel
function (DL::DoubleLayerKernel{T,S})(target,source)::T where {T,S<:Maxwell}
    k  = parameters(pde(DL))
    x  = coords(target)
    y  = coords(source)
    ny = normal(source)
    rvec = x - y
    r      = norm(rvec)
    r == 0 && return zero(T)
    g      = exp(im*k*r)/(4π*r)
    gp     = im*k*g - g/r
    rcross = cross_product_matrix(rvec)
    ncross = cross_product_matrix(ny)
    return ncross * rcross * gp/r |> transpose
end

"""
    _curl_y_green_tensor_maxwell(x, y, k)

Returns `∇ʸ × G(x, y)` where `G` is the Green tensor for Maxwell's equations
with wavenumber `k`.
"""
function _curl_y_green_tensor_maxwell(x, y, k)
    rvec = x - y
    r    = norm(rvec)
    g    = exp(im*k*r)/(4π*r)
    gp   = im*k*g - g/r
    rcross = cross_product_matrix(rvec)
    curl_G = -gp/r*rcross
    return curl_G
end

"""
    _curl_x_green_tensor_maxwell(x, y, k)

Returns `∇ˣ × G(x, y)` where `G` is the Green tensor for Maxwell's equations
with wavenumber `k`.
"""
function _curl_x_green_tensor_maxwell(x, y, k)
    rvec = x - y
    r = norm(rvec)
    g    = exp(im*k*r)/(4π*r)
    gp  = im*k*g - g/r
    rcross = cross_product_matrix(rvec)
    curl_G = gp/r*rcross
    return curl_G
end
