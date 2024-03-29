"""
    struct IntegralOperator{T,K,S,V} <: AbstractMatrix{T}

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
    T = return_type(k,X,Y)
    assert_concrete_type(T)
    IntegralOperator{T}(k,X,Y)
end

Base.size(iop::IntegralOperator)      = (length(iop.X), length(iop.Y))

kernel_type(iop::IntegralOperator) = kernel_type(iop.kernel)

kernel(iop::IntegralOperator) = iop.kernel

target_surface(iop::IntegralOperator) = iop.X
source_surface(iop::IntegralOperator) = iop.Y

Base.getindex(iop::IntegralOperator,i::Integer,j::Integer)  = getindex(kernel_type(iop),iop,i,j)

function Base.getindex(::SingleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,w = qnodes(iop.X)[i], qnodes(iop.Y)[j], qweights(iop.Y)[j]
    return iop.kernel(x,y)*w
end
function Base.getindex(::DoubleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,ny,w = qnodes(iop.X)[i], qnodes(iop.Y)[j], qnormals(iop.Y)[j], qweights(iop.Y)[j]
    return iop.kernel(x,y,ny)*w
end
function Base.getindex(::AdjointDoubleLayer,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,nx,w = qnodes(iop.X)[i], qnodes(iop.Y)[j], qnormals(iop.X)[i], qweights(iop.Y)[j]
    return iop.kernel(x,y,nx)*w
end
function Base.getindex(::HyperSingular,iop::IntegralOperator,i::Integer,j::Integer)
    x,y,nx,ny,w = qnodes(iop.X)[i], qnodes(iop.Y)[j], qnormals(iop.X)[i], qnormals(iop.Y)[j], qweights(iop.Y)[j]
    return iop.kernel(x,y,nx,ny)*w
end

combined_field_coefficients(iop::IntegralOperator) = combined_field_coefficients(iop.kernel)

# convenience constructors
SingleLayerOperator(op::AbstractPDE,X,Y=X)        = IntegralOperator(SingleLayerKernel(op),X,Y)
DoubleLayerOperator(op::AbstractPDE,X,Y=X)        = IntegralOperator(DoubleLayerKernel(op),X,Y)
AdjointDoubleLayerOperator(op::AbstractPDE,X,Y=X) = IntegralOperator(AdjointDoubleLayerKernel(op),X,Y)
HyperSingularOperator(op::AbstractPDE,X,Y=X)      = IntegralOperator(HyperSingularKernel(op),X,Y)

ambient_dimension(iop::IntegralOperator) = ambient_dimension(iop.kernel)
 

function isinside(x::SVector,mesh::NystromMesh)
    N   = ambient_dimension(mesh)     
    pde = Laplace(dim=N)
    K   = DoubleLayerKernel(pde)
    y   = qnodes(mesh)
    ν   = qnormals(mesh)   
    w   = qweights(mesh)  
    u   = sum(zip(y,ν,w)) do (y,ν,w)
         K(x,y,ν)*w
    end  
    u + 0.5 < 0 
end
isinside(x::Tuple,mesh::NystromMesh) = isinside(SVector(x),mesh)

qnodes(vec::Array{<:SVector}) = vec