# Single Layer
function (SL::SingleLayerKernel{T,S})(target,source)::T  where {T,S<:Elastostatic}
    N   = ambient_dimension(pde(SL))
    μ,λ = parameters(pde(SL))
    ν = λ/(2*(μ+λ))
    x = coords(target)
    y = coords(source)
    r = x .- y
    d = norm(r)
    d == 0 && return zero(T)
    RRT = r*transpose(r) # r ⊗ rᵗ
    if N==2
        ID = SMatrix{2,2,Float64,4}(1,0,0,1)
        return 1/(8π*μ*(1-ν))*(-(3-4*ν)*log(d)*ID + RRT/d^2)
        # return (λ + 3μ)/(4*π*(N-1)*μ*(λ+2μ))* (-log(d)*ID + (λ+μ)/(λ+3μ)*RRT/d^2)
    elseif N==3
        ID = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        return 1/(16π*μ*(1-ν)*d)*((3-4*ν)*ID + RRT/d^2)
    end
end

# Double Layer Kernel
function (DL::DoubleLayerKernel{T,S})(target,source)::T where {T,S<:Elastostatic}
    N = ambient_dimension(pde(DL))
    μ,λ = parameters(pde(DL))
    ν = λ/(2*(μ+λ))
    x = coords(target)
    y = coords(source)
    ny = normal(source)
    ν = λ/(2*(μ+λ))
    r = x .- y
    d = norm(r)
    d == 0 && return zero(T)
    RRT = r*transpose(r) # r ⊗ rᵗ
    drdn = -dot(r,ny)/d
    if N==2
        ID = SMatrix{2,2,Float64,4}(1,0,0,1)
        return -1/(4π*(1-ν)*d)*(drdn*((1-2ν)*ID + 2*RRT/d^2) + (1-2ν)/d*(r*transpose(ny) - ny*transpose(r)))
    elseif N==3
        ID = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        return -1/(8π*(1-ν)*d^2)*(drdn * ((1-2*ν)*ID + 3*RRT/d^2) + (1-2*ν)/d*(r*transpose(ny) - ny*transpose(r)))
    end
end

# Adjoint Double Layer Kernel
function (ADL::AdjointDoubleLayerKernel{T,S})(target,source)::T where {T,S<:Elastostatic}
    # DL = DoubleLayerKernel{T}(pde(ADL))
    # return -DL(x,y,nx) |> transpose
    N   = ambient_dimension(pde(ADL))
    μ,λ = parameters(pde(ADL))
    ν = λ/(2*(μ+λ))
    x = coords(target)
    nx = normal(target)
    y = coords(source)
    ν = λ/(2*(μ+λ))
    r = x .- y
    d = norm(r)
    d == 0 && return zero(T)
    RRT = r*transpose(r) # r ⊗ rᵗ
    drdn = -dot(r,nx)/d
    if N==2
        ID = SMatrix{2,2,Float64,4}(1,0,0,1)
        return -1/(4π*(1-ν)*d)*(drdn*((1-2ν)*ID + 2*RRT/d^2) + (1-2ν)/d*(r*transpose(nx) - nx*transpose(r)))
    elseif N==3
        ID = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        out =  -1/(8π*(1-ν)*d^2)*(drdn * ((1-2*ν)*ID + 3*RRT/d^2) + (1-2*ν)/d*(r*transpose(nx) - nx*transpose(r)))
        return -transpose(out)
    end
end

# Hypersingular kernel
function (HS::HyperSingularKernel{T,S})(target,source)::T where {T,S<:Elastostatic}
    N = ambient_dimension(pde(HS))
    μ,λ = parameters(pde(HS))
    ν = λ/(2*(μ+λ))
    x = coords(target)
    nx = normal(target)
    y = coords(source)
    ny = normal(source)
    r = x .- y
    d = norm(r)
    d == 0 && return zero(T)
    RRT   = r*transpose(r) # r ⊗ rᵗ
    drdn  = dot(r,ny)/d
    if N==2
        ID = SMatrix{2,2,Float64,4}(1,0,0,1)
        return μ/(2π*(1-ν)*d^2)* (2*drdn/d*( (1-2ν)*nx*transpose(r) + ν*(dot(r,nx)*ID + r*transpose(nx)) - 4*dot(r,nx)*RRT/d^2 ) +
                                  2*ν/d^2*(dot(r,nx)*ny*transpose(r) + dot(nx,ny)*RRT) +
                                  (1-2*ν)*(2/d^2*dot(r,nx)*r*transpose(ny) + dot(nx,ny)*ID + ny*transpose(nx)) -
                                  (1-4ν)*nx*transpose(ny)
                                  )
    elseif N==3
        ID = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        return μ/(4π*(1-ν)*d^3)* (3*drdn/d*( (1-2ν)*nx*transpose(r) + ν*(dot(r,nx)*ID + r*transpose(nx)) - 5*dot(r,nx)*RRT/d^2 ) +
                                  3*ν/d^2*(dot(r,nx)*ny*transpose(r) + dot(nx,ny)*RRT) +
                                  (1-2*ν)*(3/d^2*dot(r,nx)*r*transpose(ny) + dot(nx,ny)*ID + ny*transpose(nx)) -
                                  (1-4ν)*nx*transpose(ny)
                                  )
    end
end
