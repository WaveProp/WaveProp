"""
    abstract type AbstractSingularityHandler{R}
    
Used for handling localized integrand singularities in `R`.
"""
abstract type AbstractSingularityHandler{T} end

"""
    struct IMT{A,P} <: AbstractSingularityHandler{ReferenceLine}
    
One-dimensional change of variables mapping `[0,1] -> [0,1]` with the property that 
all derivatives vanish at the point `x=0`.

See [Davis and Rabinowitz](https://www.elsevier.com/books/methods-of-numerical-integration/davis/978-0-12-206360-2).
"""
struct IMT{A,P} <: AbstractSingularityHandler{ReferenceLine}
end
IMT(;a=1,p=1) = IMT{a,p}()

domain(::IMT) = ReferenceLine()
range(::IMT)  = ReferenceLine()

function (f::IMT{A,P})(x) where {A,P}
    exp(A*(1-1/x[1]^P))
end

derivative(f::IMT{A,P},x) where {A,P} = f(x) * A*P*1 / x[1]^(P+1)

jacobian(f::IMT,x) = derivative(f,x) |> SMatrix{1,1}

"""
    struct Kress{P} <: AbstractSingularityHandler{ReferenceLine}
    
Change of variables mapping `[0,1]` to `[0,1]` with the property that the first
`P` derivatives of the transformation vanish at `x=0`.
"""
struct Kress{P} <: AbstractSingularityHandler{ReferenceLine}
end
Kress(;order=5) = Kress{order}()

domain(k::Kress) = ReferenceLine()
range(k::Kress)  = ReferenceLine()

# NOTE: fastmath is needed here to allow for various algebraic simplifications
# which are not exact in floating arithmetic. Maybe reorder the operations *by
# hand* to avoid having to use fastmath? In any case, benchmark first.
@fastmath function (f::Kress{P})(y) where {P}
    x = y[1]
    v = (x) -> (1/P - 1/2)*((1-x))^3 + 1/P*((x-1)) + 1/2
    return 2v(x)^P / (v(x)^P + v(2-x)^P)
end

@fastmath function derivative(f::Kress{P},y) where {P}
    x = y[1]
    v = (x) -> (1/P - 1/2)*((1-x))^3 + 1/P*((x-1)) + 1/2
    vp = (x) -> -3*(1/P - 1/2)*((1-x))^2 + 1/P
    return 2 * (P*v(x)^(P-1)*vp(x) * (v(x)^P + v(2-x)^P) - (P*v(x)^(P-1)*vp(x) - P*v(2-x)^(P-1)*vp(2-x) ) * v(x)^P ) /
        (v(x)^P + v(2-x)^P)^2
end

jacobian(f::Kress,x) = derivative(f,x) |> SMatrix{1,1}


"""
    struct KressP{P} <: SingularityHandler{ReferenceLine}
    
Like [`Kress`](@ref), this change of variables maps the interval `[0,1]` onto
itself, but the first `P` derivatives of the transformation vanish at **both**
endpoints. 

This change of variables can be used to *periodize* integrals in the following
sense. Suppose we wish to compute the integral of `f(x)` from `0` to `1` where
`f is not a `1-periodic function. If `ϕ` is an object of type `KressP`, then
using it as a change of variables in the integration yields a similar integral
from `0` to `1` (the interval `0≤0≤1` is mappend onto itself), but with
integrand given by `g(x) = f(ϕ(x))ϕ'(x)`. Since `ϕ'` vanishes (together with `P`
of its derivatives), the function `g(x)` is now periodic (up to derivatives of
order up to `P`) at the endpoints. Thus quadrature rules designed for periodic
functions like the [`TrapezoidalP`](@ref) can be used to obtain high order
convergence of `g`, which in turn yields a modified quadrature rule when viewed
as a quadrature rule for `f`. 
"""
struct KressP{P} <: AbstractSingularityHandler{ReferenceLine}
end
KressP(;order=5) = KressP{order}()

domain(k::KressP) = ReferenceLine()
range(k::KressP)  = ReferenceLine()

@fastmath function (f::KressP{P})(y) where {P}
    x = y[1]
    v = (x) -> (1/P - 1/2)*((1-2x))^3 + 1/P*((2x-1)) + 1/2
    return v(x)^P / (v(x)^P + v(1-x)^P)
end

@fastmath function derivative(f::KressP{P},y) where {P}
    x = y[1]
    v =  (x) -> (1/P - 1/2)*((1-2x))^3 + 1/P*((2x-1)) + 1/2    
    vp = (x) -> -6*(1/P - 1/2)*((1-2x))^2 + 2/P
    return (P*v(x)^(P-1)*vp(x) * (v(x)^P + v(1-x)^P) - (P*v(x)^(P-1)*vp(x) - P*v(1-x)^(P-1)*vp(1-x) ) * v(x)^P ) /
        (v(x)^P + v(1-x)^P)^2
end

jacobian(f::KressP,x) = derivative(f,x) |> SMatrix{1,1}



"""
    struct Duffy <: AbstractSingularityHandler{RefereceTriangle}

Change of variables mapping the `ReferenceSquare` to the `RefereceTriangle` with
the property that the jacobian vanishes at the `(1,0)` vertex of the triangle.

Useful for integrating functions with a singularity on the `(1,0)` edge of the
reference triangle.
"""
struct Duffy <: AbstractSingularityHandler{ReferenceTriangle} end

domain(::Duffy) = ReferenceSquare()
range(::Duffy)  = ReferenceTriangle()

function (::Duffy)(u)
    SVector(u[1][1],(1-u[1][1])*u[2][1])
end    

function jacobian(::Duffy,u)
    SMatrix{2,2,Float64}(1,0,-u[2][1],(1-u[1][1]))
end    

# TODO: generalize to `N` dimensions
"""
    struct TensorProductSingularityHandler{S} <:
    AbstractSingularityHandler{ReferenceSquare}
        
A tensor product of two one-dimensional `SingularityHandler`s for performing
integration over the `ReferenceSquare`.
"""
struct TensorProductSingularityHandler{S} <: AbstractSingularityHandler{ReferenceSquare} 
    shandler::S
end

domain(::TensorProductSingularityHandler) = ReferenceSquare()

function TensorProductSingularityHandler(q...)
    TensorProductSingularityHandler(q)
end    

function (f::TensorProductSingularityHandler)(x)
    shandler = f.shandler    
    @assert length(shandler) == length(x)
    N = length(shandler)
    svector(i->shandler[i](x[i]),N)
end    

function jacobian(f::TensorProductSingularityHandler,x)
    shandler = f.shandler    
    @assert length(shandler) == length(x)
    N = length(shandler)
    if N == 2
        jx = derivative(shandler[1],x[1])
        jy = derivative(shandler[2],x[2])
        return SMatrix{2,2}(jx,0,0,jy)
    else 
        notimplemented()
    end
end    

# # overload quadgk so that a singularity handler can be passed as first argument
# function quadgk(phi::SingularityHandler{ReferenceLine},f,args...;xs,kwargs...)
#     a,b = args[1], args[end]
#     quadgk(args...;kwargs...) do x
#         xi,dxi = phi(x,a,b,xs), derivative(phi,x,a,b,xs)
#         f(xi)*dxi
#     end
# end

# """
#     Window{A,B,S}

# Change of variables mapping `[0,1]` to `[0,1]` with the following properties:
# - smooth (infinitely differentiable) on `[0,1]`
# - all derivatives vanish on both `x=0` and `x=1`
# - derivative of the transoformation is exactly 1 between `A` and `B`
# """
# struct Window{A,B,S} <: AbstractSingularityHandler{ReferenceLine}
# end
# Window() = Window{1,1,5}()

# domain(::Window) = ReferenceLine()
# range(::Window)  = ReferenceLine()

# normalization(w::Window) = quadgk(t->_derivative(w,t),0,1,rtol=1e-16)[1]

# function (f::Window{A,B,S})(x) where {A,B,S}
#     x == 0 && return 0.0
#     I,_ = quadgk(t->derivative(f,t),0,x)  
#     return I  
#     # FIXME: the current way of computing the values of the `Window` change of
#     # variables is very innefficient. Could use a e.g. `cumsum`. 
# end    

# function derivative(f::Window,x)
#     x==0 && (return 0.0) 
#     x==1 && (return 0.0)
#     _derivative(f,x) / normalization(f)
# end    

# function _derivative(f::Window{A,B,S},x) where {A,B,S}
#     if  A ≤ x ≤ B
#         return 1
#     elseif 0 ≤ x < A
#         u = (A - x) / (A)
#         return exp(S*exp(-1/u)/(u-1))
#     elseif B < x ≤ 1
#         u = (x - B) / (1 - B)
#         return exp(S*exp(-1/u)/(u-1))
#     else
#         return 0
#     end
# end    

# jacobian(f::Window,x) = derivative(f,x) |> SMatrix{1,1}