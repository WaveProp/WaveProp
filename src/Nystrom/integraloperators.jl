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

Base.size(iop::IntegralOperator)      = (length(qnodes(iop.X)), length(qnodes(iop.Y)))

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

function SingularIntegration.singular_weights(iop::IntegralOperator,qstd::AbstractQuadratureRule,q::SingularQuadratureRule)
    X,Y = iop.X, iop.Y    
    T = eltype(iop)    
    K = iop.kernel
    dict_near = near_interaction_list(X,Y;dim=ambient_dimension(Y)-1,atol=0.1)
    Is    = Int[]
    Js    = Int[]
    Vs    = T[]
    ui    = qstd()[1]
    for (E,list_near) in dict_near # for each element type
        el2qnodes = Y.el2qnodes[E]
        num_qnodes, num_els   = size(el2qnodes)
        for (n,el) in enumerate(ElementIterator{E}(Y)) # for each element of type E in the mesh Y
            j_glob                = el2qnodes[:,n]
            for (i,jloc) in list_near[n] # for each near observation point
                x   = X.qnodes[i] # observation point
                vs  = ui[jloc]    # regularization center in parameter space
                if kernel_type(K) == SingleLayer()
                    k    = (v) -> K(x,el(v))*measure(el,v)
                elseif kernel_type(K) == DoubleLayer()
                    k    = (v) -> K(x,el(v),normal(el,v)) * measure(el,v)
                end    
                w  = singular_weights(k,ui,q,vs)
                append!(Is,fill(i,num_qnodes))
                append!(Js,j_glob)
                w  = w - iop[i,j_glob]
                append!(Vs,w)
            end            
        end
    end  
    return Sp = sparse(Is,Js,Vs,size(iop)...)      
end    