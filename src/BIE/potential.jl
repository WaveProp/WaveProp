struct IntegralPotential{T<:AbstractKernel,S}
    kernel::T
    surface::S
end

function IntegralPotential(op::AbstractPDE,surf,ktype::Symbol)
    if ktype == :singlelayer
        k = SingleLayerKernel(op)
    elseif ktype == :doublelayer
        k = DoubleLayerKernel(op)
    else
        error("unknown kernel type: valid options are :singlelayer or :doublelayer")
    end
    return IntegralPotential(k,surf)
end

kernel_type(pot::IntegralPotential) = kernel_type(pot.kernel)

(pot::IntegralPotential)(σ::AbstractVector,x) = pot(kernel_type(pot),σ,x)
function (pot::IntegralPotential)(::SingleLayer,σ::AbstractVector,x)
    f = pot.kernel
    Γ = pot.surface
    iter = zip(qnodes(Γ),qweights(Γ),σ)
    out = mapreduce(+,iter) do (y,w,σ)
        f(x,y)*σ*w
    end
    return out
end
function (pot::IntegralPotential)(::DoubleLayer,σ::AbstractVector,x)
    f = pot.kernel
    Y = pot.surface
    iter = zip(qnodes(Y),qweights(Y),qnormals(Y),σ)
    out  = mapreduce(+,iter) do (y,w,ny,σ)
        f(x,y,ny)*σ*w
    end
    return out
end

Base.getindex(pot::IntegralPotential,σ::AbstractVector) = (x) -> pot(σ,x)

SingleLayerPotential(op::AbstractPDE,surf) = IntegralPotential(SingleLayerKernel(op),surf)
DoubleLayerPotential(op::AbstractPDE,surf) = IntegralPotential(DoubleLayerKernel(op),surf)
