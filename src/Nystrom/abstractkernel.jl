"""
    abstract type AbstractKernel{T}

A kernel functions `K` with the signature `K(target,source)::T`.
"""
abstract type AbstractKernel{T} end

return_type(K::AbstractKernel{T}) where {T} = T

ambient_dimension(K::AbstractKernel) = ambient_dimension(pde(K))

pde(k::AbstractKernel) = k.pde

"""
    struct GenericKernel{T,F} <: AbstractKernel{T}

An [`AbstractKernel`](@ref) with `kernel` of type `F`.
"""
struct GenericKernel{T,F} <: AbstractKernel{T}
    kernel::F
end

"""
    struct SingleLayerKernel{T,Op} <: AbstractKernel{T}

Given an operator `Op`, construct its free-space single-layer kernel (i.e. the
fundamental solution).
"""
struct SingleLayerKernel{T,Op} <: AbstractKernel{T}
    pde::Op
end
SingleLayerKernel{T}(op) where {T} = SingleLayerKernel{T,typeof(op)}(op)
SingleLayerKernel(op)              = SingleLayerKernel{default_kernel_eltype(op)}(op)

"""
    struct DoubleLayerKernel{T,Op} <: AbstractKernel{T}

Given an operator `Op`, construct its free-space double-layer kernel. This
corresponds to the `γ₁` trace of the [`SingleLayerKernel`](@ref). For operators
such as [`Laplace`](@ref) or [`Helmholtz`](@ref), this is simply the normal
derivative of the fundamental solution respect to the source variable.
"""
struct DoubleLayerKernel{T,Op} <: AbstractKernel{T}
    pde::Op
end
DoubleLayerKernel{T}(op) where {T}     = DoubleLayerKernel{T,typeof(op)}(op)
DoubleLayerKernel(op)                  = DoubleLayerKernel{default_kernel_eltype(op)}(op)

"""
    struct AdjointDoubleLayerKernel{T,Op} <: AbstractKernel{T}

Given an operator `Op`, construct its free-space adjoint double-layer kernel. This
corresponds to the `transpose(γ₁,ₓ[G])`, where `G` is the
[`SingleLayerKernel`](@ref). For operators such as [`Laplace`](@ref) or
[`Helmholtz`](@ref), this is simply the normal derivative of the fundamental
solution respect to the target variable.
"""
struct AdjointDoubleLayerKernel{T,Op} <: AbstractKernel{T}
    pde::Op
end
AdjointDoubleLayerKernel{T}(op) where {T} = AdjointDoubleLayerKernel{T,typeof(op)}(op)
AdjointDoubleLayerKernel(op)              = AdjointDoubleLayerKernel{default_kernel_eltype(op)}(op)

"""
    struct HyperSingularKernel{T,Op} <: AbstractKernel{T}

Given an operator `Op`, construct its free-space hypersingular kernel. This
corresponds to the `transpose(γ₁,ₓγ₁[G])`, where `G` is the
[`SingleLayerKernel`](@ref). For operators such as [`Laplace`](@ref) or
[`Helmholtz`](@ref), this is simply the normal derivative of the fundamental
solution respect to the target variable of the `DoubleLayerKernel`.
"""
struct HyperSingularKernel{T,Op} <: AbstractKernel{T}
    pde::Op
end
HyperSingularKernel{T}(op) where {T} = HyperSingularKernel{T,typeof(op)}(op)
HyperSingularKernel(op)              = HyperSingularKernel{default_kernel_eltype(op)}(op)

# a trait for the kernel type
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
