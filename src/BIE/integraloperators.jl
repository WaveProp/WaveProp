"""
    IntegralOperator{T,K,S,V} <: AbstractMatrix{T}

Representation of an integral operator which takes a density ϕ defined on `Y`
and integrates it with `kernel` for all elements `x ∈ X`.
"""
struct IntegralOperator{T,K,S,V} <: AbstractMatrix{T}
    kernel::K
    X::S
    Y::V
end
IntegralOperator{T}(k,X,Y) where {T} = IntegralOperator{T,typeof(k),typeof(X),typeof(Y)}(k,X,Y)
function IntegralOperator(k,X,Y=X)
    T = return_type(k)
    IntegralOperator{T}(k,X,Y)
end

Base.size(iop::IntegralOperator)      = (length(nodes(iop.X)), length(nodes(iop.Y)))

kernel_type(iop::IntegralOperator) = kernel_type(iop.kernel)

Base.getindex(iop::IntegralOperator,i::Integer,j::Integer)  = getindex(kernel_type(iop),iop,i,j)

function Base.getindex(::SingleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,w = nodes(iop.X)[i], nodes(iop.Y)[j], weights(iop.Y)[j]
    return iop.kernel(x,y)*w
end
function Base.getindex(::DoubleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,ny,w = nodes(iop.X)[i], nodes(iop.Y)[j], normals(iop.Y)[j], weights(iop.Y)[j]
    return iop.kernel(x,y,ny)*w
end
function Base.getindex(::AdjointDoubleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,nx,w = nodes(iop.X)[i], nodes(iop.Y)[j], normals(iop.X)[i], weights(iop.Y)[j]
    return iop.kernel(x,y,nx)*w
end
function Base.getindex(::HyperSingular,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,nx,ny,w = nodes(iop.X)[i], nodes(iop.Y)[j], normals(iop.X)[i], normals(iop.Y)[j], weights(iop.Y)[j]
    return iop.kernel(x,y,nx,ny)*w
end

# convenience constructors
SingleLayerOperator(op::AbstractPDE,X,Y=X)        = IntegralOperator(SingleLayerKernel(op),X,Y)
DoubleLayerOperator(op::AbstractPDE,X,Y=X)        = IntegralOperator(DoubleLayerKernel(op),X,Y)
AdjointDoubleLayerOperator(op::AbstractPDE,X,Y=X) = IntegralOperator(AdjointDoubleLayerKernel(op),X,Y)
HyperSingularOperator(op::AbstractPDE,X,Y=X)      = IntegralOperator(HyperSingularKernel(op),X,Y)

ambient_dimension(iop::IntegralOperator) = ambient_dimension(iop.kernel)
