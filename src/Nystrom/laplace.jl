# Single Layer
function (SL::SingleLayerKernel{T,Laplace{N}})(target,source)::T  where {N,T}
    x = coords(target)
    y = coords(source)
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N==2
        return -1/(2π)*log(d)
    elseif N==3
        return 1/(4π)/d
    else
        notimplemented()
    end
end

# Double Layer Kernel
function (DL::DoubleLayerKernel{T,Laplace{N}})(target,source)::T where {N,T}
    x = coords(target)
    y = coords(source)
    ny = normal(source)
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N == 2
        return 1/(2π)/(d^2) .* dot(r,ny)
    elseif N==3
        return 1/(4π)/(d^3) .* dot(r,ny)
    else
        notimplemented()
    end
end

# Adjoint double Layer Kernel
function (ADL::AdjointDoubleLayerKernel{T,Laplace{N}})(target,source)::T where {N,T}
    x = coords(target)
    y = coords(source)
    nx = normal(target)
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N==2
        return -1/(2π)/(d^2) .* dot(r,nx)
    elseif N==3
        return -1/(4π)/(d^3) .* dot(r,nx)
    end
end

# Hypersingular kernel
function (HS::HyperSingularKernel{T,Laplace{N}})(target,source)::T where {N,T}
    x = coords(target)
    y = coords(source)
    nx = normal(target)
    ny = normal(source)
    r = x - y
    d = norm(r)
    d == 0 && (return zero(T))
    if N==2
        ID = Mat{2,2,Float64,4}(1,0,0,1)
        RRT = r*transpose(r) # r ⊗ rᵗ
        return 1/(2π)/(d^2) * transpose(nx)*(( ID -2*RRT/d^2  )*ny)
    elseif N==3
        ID = Mat{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        RRT = r*transpose(r) # r ⊗ rᵗ
        return 1/(4π)/(d^3) * transpose(nx)*(( ID -3*RRT/d^2  )*ny)
    end
end
