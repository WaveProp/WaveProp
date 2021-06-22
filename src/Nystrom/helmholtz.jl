function (SL::SingleLayerKernel{T,S})(target, source)::T where {T,S <: Helmholtz}
    x = coords(target)
    y = coords(source)
    N = ambient_dimension(pde(SL))
    k = parameters(pde(SL))
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N == 2
        return im / 4 * hankelh1(0, k * d)
    elseif N == 3
        return 1 / (4π) / d * exp(im * k * d)
    end
end

# Double Layer Kernel
function (DL::DoubleLayerKernel{T,S})(target, source)::T where {T,S <: Helmholtz}
    x, y, ny = coords(target), coords(source), normal(source)
    N = ambient_dimension(pde(DL))
    k = parameters(pde(DL))
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N == 2
        return im * k / 4 / d * hankelh1(1, k * d) .* dot(r, ny)
    elseif N == 3
        return 1 / (4π) / d^2 * exp(im * k * d) * ( -im * k + 1 / d ) * dot(r, ny)
    end
end

# Adjoint double Layer Kernel
function (ADL::AdjointDoubleLayerKernel{T,S})(target, source)::T where {T,S <: Helmholtz}
    x, y, nx = coords(target), coords(source), normal(target)
    N = ambient_dimension(pde(ADL))
    k = parameters(pde(ADL))
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N == 2
        return -im * k / 4 / d * hankelh1(1, k * d) .* dot(r, nx)
    elseif N == 3
        return -1 / (4π) / d^2 * exp(im * k * d) * ( -im * k + 1 / d ) * dot(r, nx)
    end
end

# Hypersingular kernel
function (HS::HyperSingularKernel{T,S})(target, source)::T where {T,S <: Helmholtz}
    x, y, nx, ny = coords(target), coords(source), normal(target), normal(source)
    N = ambient_dimension(pde(HS))
    k = parameters(pde(HS))
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N == 2
        RRT = r * transpose(r) # r ⊗ rᵗ
        # TODO: rewrite the operation below in a more clear/efficient way
        return transpose(nx) * ((-im * k^2 / 4 / d^2 * hankelh1(2, k * d) * RRT + im * k / 4 / d * hankelh1(1, k * d) * I) * ny)
    elseif N == 3
        RRT   = r * transpose(r) # r ⊗ rᵗ
        term1 = 1 / (4π) / d^2 * exp(im * k * d) * ( -im * k + 1 / d ) * I
        term2 = RRT / d * exp(im * k * d) / (4 * π * d^4) * (3 * (d * im * k - 1) + d^2 * k^2)
        return  transpose(nx) * (term1 + term2) * ny
    end
end
