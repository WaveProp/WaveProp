function (SL::SingleLayerKernel{T,S})(x,y)::T  where {T,S<:Maxwell}
    N = ambient_dimension(S)    
    x==y && return zero(T)
    k = SL.op.k
    r = x - y
    d = norm(r)
    if N==2
        return error("Maxwell operator not implemented in 2d")
    elseif N==3
        g   = 1/(4π)/d * exp(im*k*d)
        gp  = im*k*g - g/d
        gpp = im*k*gp - gp/d + g/d^2
        ID    = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        RRT   = r*transpose(r) # r ⊗ rᵗ
        return  g*ID + 1/k^2*(gp/d*ID + (gpp/d^2 - gp/d^3)*RRT)
    end
end

# Double Layer Kernel
# n × ∇ × G = γ₁ G
function (DL::DoubleLayerKernel{T,S})(x,y,ny)::T where {T,S<:Maxwell}
    N = ambient_dimension(S)    
    x==y && return zero(T)
    k = DL.op.k
    r = x - y
    d = norm(r)
    g   = 1/(4π)/d * exp(im*k*d)
    gp  = im*k*g - g/d
    if N==2
        return error("Maxwell operator not yet defined in 2d")
    elseif N==3
        ID    = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
        ncross = transpose(SMatrix{3,3,Float64,9}(0,-ny[3],ny[2],
                                         ny[3],0,-ny[1],
                                         -ny[2],ny[1],0))
        rcross = transpose(SMatrix{3,3,Float64,9}(0,-r[3],r[2],
                                         r[3],0,-r[1],
                                         -r[2],r[1],0))
        # return -gp/d*ncross*rcross
        return gp/d*rcross*ncross
    end
end
