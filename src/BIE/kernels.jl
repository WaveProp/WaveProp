abstract type AbstractKernel{T} end

return_type(K::AbstractKernel{T}) where {T} = T
ambient_dimension(K::AbstractKernel)              = ambient_dimension(K.op)

struct SingleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
SingleLayerKernel{T}(op) where {T} = SingleLayerKernel{T,typeof(op)}(op)
SingleLayerKernel(op)              = SingleLayerKernel{default_kernel_eltype(op)}(op)

struct DoubleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
DoubleLayerKernel{T}(op) where {T} = DoubleLayerKernel{T,typeof(op)}(op)
DoubleLayerKernel(op)              = DoubleLayerKernel{default_kernel_eltype(op)}(op)

struct AdjointDoubleLayerKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
AdjointDoubleLayerKernel{T}(op) where {T} = AdjointDoubleLayerKernel{T,typeof(op)}(op)
AdjointDoubleLayerKernel(op)              = AdjointDoubleLayerKernel{default_kernel_eltype(op)}(op)

struct HyperSingularKernel{T,Op} <: AbstractKernel{T}
    op::Op
end
HyperSingularKernel{T}(op) where {T} = HyperSingularKernel{T,typeof(op)}(op)
HyperSingularKernel(op)              = HyperSingularKernel{default_kernel_eltype(op)}(op)

kernel_type(::SingleLayerKernel)        = SingleLayer()
kernel_type(::DoubleLayerKernel)        = DoubleLayer()
kernel_type(::AdjointDoubleLayerKernel) = AdjointDoubleLayer()
kernel_type(::HyperSingularKernel)      = HyperSingular()