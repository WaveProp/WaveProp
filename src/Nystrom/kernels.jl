abstract type AbstractKernel{T} end

return_type(K::AbstractKernel{T},args...) where {T} = T
ambient_dimension(K::AbstractKernel) = ambient_dimension(K.op)

struct SingleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
SingleLayerKernel{T}(op) where {T} = SingleLayerKernel{T,typeof(op)}(op)
SingleLayerKernel(op)              = SingleLayerKernel{default_kernel_eltype(op)}(op)
(op::SingleLayerKernel)(x,y)       = op(SVector(x),SVector(y))

struct DoubleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
DoubleLayerKernel{T}(op) where {T} = DoubleLayerKernel{T,typeof(op)}(op)
DoubleLayerKernel(op)              = DoubleLayerKernel{default_kernel_eltype(op)}(op)
(op::DoubleLayerKernel)(x,y,ny)       = op(SVector(x),SVector(y),SVector(ny))

struct AdjointDoubleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
AdjointDoubleLayerKernel{T}(op) where {T} = AdjointDoubleLayerKernel{T,typeof(op)}(op)
AdjointDoubleLayerKernel(op)              = AdjointDoubleLayerKernel{default_kernel_eltype(op)}(op)
(op::AdjointDoubleLayerKernel)(x,y,nx)       = op(SVector(x),SVector(y),SVector(nx))

struct HyperSingularKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
HyperSingularKernel{T}(op) where {T} = HyperSingularKernel{T,typeof(op)}(op)
HyperSingularKernel(op)              = HyperSingularKernel{default_kernel_eltype(op)}(op)
(op::HyperSingularKernel)(x,y,nx,ny) = op(SVector(x),SVector(y),SVector(nx),SVector(ny))

# kernel trait do distinguish various integral operators and potentials. This
# helps in dispatch and in known the expected kernel signature (e.g. does it
# require the normal at point x? What about at the ponit y? etc)
struct SingleLayer end
struct DoubleLayer end
struct AdjointDoubleLayer end
struct HyperSingular end

kernel_type(::SingleLayerKernel)        = SingleLayer()
kernel_type(::DoubleLayerKernel)        = DoubleLayer()
kernel_type(::AdjointDoubleLayerKernel) = AdjointDoubleLayer()
kernel_type(::HyperSingularKernel)      = HyperSingular()

combined_field_coefficients(::SingleLayerKernel)        = (0,-1)
combined_field_coefficients(::DoubleLayerKernel)        = (1,0)
combined_field_coefficients(::AdjointDoubleLayerKernel) = (0,-1)
combined_field_coefficients(::HyperSingularKernel)      = (1,0)

# most of the types, you want to specifically define this method for a given
# kernel K. In the generic case, however, we rely on promote_op to infer the
# return type given the trait `kernel_type`
function return_type(K,X,Y)
    x_type = eltype(qnodes(X))
    y_type = eltype(qnodes(Y))
    if kernel_type(K) == SingleLayer()
        return Base.promote_op(K,x_type,y_type)
    elseif kernel_type(K) == DoubleLayer()
        νy_type = eltype(qnormals(Y))
        return Base.promote_op(K,x_type,y_type,νy_type)
    elseif kernel_type(K) == AdjointDoubleLayer()
        νx_type = eltype(qnormals(X))    
        return Base.promote_op(K,x_type,y_type,νx_type)
    elseif kernel_type(K) == HyperSingular()
        νx_type = eltype(qnormals(X))
        νy_type = eltype(qnormals(Y))
        return Base.promote_op(K,x_type,y_type,νx_type,νy_type)
    else
        error("unknown `kernel_type` trait")        
    end
end    

