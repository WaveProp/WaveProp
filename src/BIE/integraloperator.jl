"""
    IntegralOperator{V,Tx,Ty,Tk} <: AbstractMatrix{V}

# Fields
- `kernel::K`
- `target::T`
- `source::S`

An integral operator with `kernel::K` mapping a density defined on `source` to a density defined on `target`.
"""
mutable struct IntegralOperator{V,K,T,S} <: AbstractMatrix{V}
    kernel::K
    target::T
    source::S
end

target(iop::IntegralOperator) = iop.target
source(iop::IntegralOperator) = iop.source

IntegralOperator{V}(kernel,X,Y=X) where {V} = IntegralOperator{V,typeof(kernel),typeof(X),typeof(Y)}(kernel,X,Y)

kerneltype(iop::IntegralOperator) = kerneltype(iop.kernel)

get_weights(iop::IntegralOperator,args...) = get_weights(source_surface(iop),args...)

# Base.size(iop::IntegralOperator) = length(target_surface(iop)),length(source_surface(iop))

# Base.getindex(iop::IntegralOperator,i::Int,j::Int) = compute_kernel_entry(iop,i,j)*getweights(iop,j)

# """
#     compute_kernel_entry(k,X,Y,i,j)

# Specify how to compute the `(i,j)` kernel entry from the function `k` and the
# meshes `X` and `Y`. This will dispatch depending on the `kerneltype` trait of `k`.
# """
# compute_kernel_entry(f,X,Y,i,j) = _compute_kernel_entry(kerneltype(f),f,X,Y,i,j)

# """
#     compute_kernel_entry(iop::IntegralOperator,i,j)

# Compute the `(i,j)` entry of the kernel part of the integral operator; that is, without the quadrature weights.

# It is useful to separate the two in cases where e.g. the kernel is symmetric. 
# """
# compute_kernel_entry(iop::IntegralOperator,i,j) = compute_kernel_entry(iop.kernel,target_surface(iop),source_surface(iop),i,j)

# function _compute_kernel_entry(::Union{SingleLayer,UnknownKernelType},f,X,Y,i,j)
#     x  = getnodes(X,i)
#     y  = getnodes(Y,j)
#     return f(x,y)
# end

# function _compute_kernel_entry(::Union{DoubleLayer},f,X,Y,i,j)
#     x  = getnodes(X,i)
#     y  = getnodes(Y,j)
#     ny = getnormals(Y,j)
#     return f(x,y,ny)
# end

# function _compute_kernel_entry(::Union{AdjointDoubleLayer},f,X,Y,i,j)
#     x  = getnodes(X,i)
#     nx = getnormals(X,j)
#     y  = getnodes(Y,j)
#     return f(x,nx,y)
# end

# function _compute_kernel_entry(::Union{HyperSingular},f,X,Y,i,j)
#     x  = getnodes(X,i)
#     nx = getnormals(X,j)
#     y  = getnodes(Y,j)
#     ny = getnormals(Y,j)
#     return f(x,nx,y,ny)
# end

# function IntegralOperator(f,X,Y=X)
#     # # try to infer return type of kernel. To be on the safe side, throw an error if inference fails
#     T = Base.promote_op(compute_kernel_entry,typeof(f),typeof(X),typeof(Y),Int,Int)
#     isconcretetype(T) || throw(ConcreteInferenceError(T))
#     return IntegralOperator{T}(f,X,Y)
# end

# function self_interaction(op::IntegralOperator)
#     X,Y = target_surface(op), source_surface(op)
#     for τ in getelements(X)
#         self_interaction_block(op,τ)
#     end
# end

# function self_interaction_block(op,τ)
#     n = length(τ)
#     B = Matrix{eltype(op)}(undef,n,n)
#     K = op.kernel
#     k = (u,v) -> K(τ(u),τ(v))*det(jacobian(τ,v))
#     for i in 1:n
#         u       = parametric_coordinate(τ,i)
#         qs      = ϕ(q,xs=u)
#         A       = lagrange_interpolation_matrix(getnodes(qs),getnodes(q))
#         B[i,:]  = transpose(A)*k.(u,nodes(qs))
#     end
# end

# """
#     compute_kernel_entry(k,X,Y,i,j)

# Specify how to compute the `(i,j)` kernel entry from the function `k` and the
# meshes `X` and `Y`.

# For the *Nystrom* boundary integral methods, this usually comes down to simply
# evaluating `k(getnode(X,i),getnode(Y,j))*getweight(Y,j)`, but other methods such
# as boundary element methods may require more information from the meshes `X` and
# `Y`.

# The implementation is responsible for defining this method.
# """
# function compute_kernel_entry
# end

# # dispatch depending on availability of fast methods
# function mul!(Y::Vector,IOp::IntegralOperator,X::Vector,a::Number,b::Number)
#     if hashmatrix(IOp)
#         mul!(Y,IOp.hmat,X,a,b)
#     elseif hasfmm(IOp)
#         mul!(Y,IOp.fmm,X,a,b)
#     elseif hasmatrix(IOp)
#         mul!(Y,IOp.matrix,X,a,b)
#     else
#         #fallback to generic on Base
#         invoke(mul!,Tuple{typeof(Y),AbstractMatrix{eltype(IOp)},typeof(X),typeof(a),typeof(b)},Y,IOp,X,a,b)
#     end
#     return Y
# end

# ## methods to complete IOp construction (not needed, but useful for e.g. acceleration)
# """
#     compute_matrix!(IOp::IntegralOperator)

# Compute a full matrix representation of `IOp` and store it in `IOp.matrix`.
# """
# function compute_matrix!(IOp::IntegralOperator)
#     IOp.matrix = Matrix(IOp)
#     return IOp
# end

# """
#     compute_correction!(op::IntegralOperator;algo)

# Compute a sparse correction to the integral operator `op`.

# The `algo` argument should be a callable object with signature
#  `algo(::IntegralOperator)::SparseCSV{V}` (e.g. a function). A default algorithm
#  [`green_correction!`](@ref) is provided which works with the currently
#  implemented kernels, but users can supply their own algorithm.

# See also: [`green_correction!`](@ref)
# """
# function compute_correction!(algo,op::IntegralOperator)
#     op.correction = algo(op)
# end

# ## convenience methods
# kerneltype(IOp::IntegralOperator) = kerneltype(IOp.kernel)
# hashmatrix(IOp::IntegralOperator) = length(IOp.hmatrix) > 0
# hasmatrix(IOp::IntegralOperator) = length(IOp.matrix) > 0
# hasfmm(IOp::IntegralOperator)    = false
# hascorrection(IOp::IntegralOperator)    = length(IOp.correction) > 0
