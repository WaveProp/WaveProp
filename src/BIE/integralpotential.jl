"""
    IntegralPotential{T,S}

A functor-like object representing the integral over `surface` of `kernel` multiplied by a given `density`.

Can be called using `pot(σ::Density,x)` to return the value of the integral potential at `x`, or using
`pot[σ]` to return an anonymous function.
"""
struct IntegralPotential{K,S}
    kernel::K
    density::Density{T,S}
    surface::S
end

# NOTE: the evaluation of IntegralPotential must first dispatch on a kernel type
# in order to know how to evaluate the kernel function. It the kerneltype trait
# is SingleLayer, then the kernels is called as k(x,y), while for the double
# layer the kernel is called as k(x,y,ny).
(pot::IntegralPotential)(σ::AbstractVector,x) = pot(kerneltype(pot),σ,x)
function (pot::IntegralPotential)(::Union{SingleLayer,UnknownKernelType},σ::AbstractVector,x)
    f = pot.kernel
    Γ = pot.surface
    iter = zip(getnodes(Γ),getweights(Γ),σ)
    out = mapreduce(+,iter) do (y,w,σ)
        f(x,y)*σ*w
    end
    return out
end
function (pot::IntegralPotential)(::DoubleLayer,σ::AbstractVector,x)
    f = pot.kernel
    Y = pot.surface
    iter = zip(getnodes(Y),getweights(Y),getnormals(Y),σ)
    out  = mapreduce(+,iter) do (y,w,ny,σ)
        f(x,y,ny)*σ*w
    end
    return out
end

Base.getindex(pot::IntegralPotential,σ::AbstractVector) = (x) -> pot(σ,x)

# function IntegralPotential(op::AbstractPDE,surf,ktype::Symbol)
#     if ktype == :singlelayer
#         k = SingleLayerKernel(op)
#     elseif ktype == :doublelayer
#         k = DoubleLayerKernel(op)
#     else
#         error("unknown kernel type: valid options are :singlelayer or :doublelayer")
#     end
#     return IntegralPotential(k,surf)
# end

# kernel_type(pot::IntegralPotential) = kernel_type(pot.kernel)

# Base.getindex(pot::IntegralPotential,σ::AbstractVector) = (x) -> pot(σ,x)

# SingleLayerPotential(op::AbstractPDE,surf) = IntegralPotential(SingleLayerKernel(op),surf)
# DoubleLayerPotential(op::AbstractPDE,surf) = IntegralPotential(DoubleLayerKernel(op),surf)

Base.getindex(pot::IntegralPotential,σ::Density) = (x) -> pot(σ,x)
