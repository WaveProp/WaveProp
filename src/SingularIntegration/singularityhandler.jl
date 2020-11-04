"""
    abstract type SingularityHandler{T}
    
Functor types used in handling singularies in the integrand when using
quadrature rules.        
"""
abstract type SingularityHandler{T} end

# generic change of variables formula
(f::SingularityHandler{ReferenceLine})(x,a,b)    = a + (b-a)*f((x-a)/(b-a))
derivative(f::SingularityHandler{ReferenceLine},x,a,b)    = derivative(f,(x-a)/(b-a))
(f::SingularityHandler{ReferenceLine})(x,a,b,xs) = x < xs ? f(x,xs,a) : f(x,xs,b)
derivative(f::SingularityHandler{ReferenceLine},x,a,b,xs)    = x < xs ? derivative(f,x,xs,a) : derivative(f,x,xs,b)

function QuadGK.quadgk(phi::SingularityHandler{ReferenceLine},f,args...;xs,kwargs...)
    a,b = args[1], args[end]
    QuadGK.quadgk(args...;kwargs...) do x
        xi,dxi = phi(x,a,b,xs), derivative(phi,x,a,b,xs)
        f(xi)*dxi
    end
end

"""
    struct IMT{A,P} <: SingularityHandler{ReferenceLine}
    
One-dimensional change of variables mapping `[0,1] -> [0,1]` with the property that 
all derivatives vanish at the point `x=0`.

# FIXME: how do you add a reference to a docstring?
See [Davis and Rabinowitz](https://www.elsevier.com/books/methods-of-numerical-integration/davis/978-0-12-206360-2)
"""
struct IMT{A,P} <: SingularityHandler{ReferenceLine}
end
IMT(;a=1,p=1) = IMT{a,p}()

domain(::IMT) = ReferenceLine()
range(::IMT)  = ReferenceLine()

function (f::IMT{A,P})(x) where {A,P}
    exp(A*(1-1/x^P))
end

derivative(f::IMT{A,P},x) where {A,P} = f(x) * A*P*1 / x^(P+1)
jacobian(f::IMT,x) = derivative(f,x) |> SMatrix{1,1}

"""
    struct Kress{P} <: SingularityHandler{ReferenceLine}
    
Change of variables mapping `[0,1]` to `[0,1]` with the property that the first
`P` derivatives of the transformation vanish at `x=0`.
"""
struct Kress{P} <: SingularityHandler{ReferenceLine}
end
Kress(;order=5) = Kress{order}()

domain(k::Kress) = ReferenceLine()
range(k::Kress)  = ReferenceLine()

# TODO: reorder operations of Kress so that fastmath is not needed
@fastmath function (f::Kress{P})(x) where {P}
    v = (x) -> (1/P - 1/2)*((1-x))^3 + 1/P*((x-1)) + 1/2
    return 2v(x)^P / (v(x)^P + v(2-x)^P)
end
jacobian(f::Kress,x) = derivative(f,x) |> SMatrix{1,1}

@fastmath function derivative(f::Kress{P},x) where {P}
    v = (x) -> (1/P - 1/2)*((1-x))^3 + 1/P*((x-1)) + 1/2
    vp = (x) -> -3*(1/P - 1/2)*((1-x))^2 + 1/P
    return 2 * (P*v(x)^(P-1)*vp(x) * (v(x)^P + v(2-x)^P) - (P*v(x)^(P-1)*vp(x) - P*v(2-x)^(P-1)*vp(2-x) ) * v(x)^P ) /
        (v(x)^P + v(2-x)^P)^2
end

struct Window{A,B,S} end

function (f::Window{A,B,S})(x) where {A,B,S}
    # TODO: Can the window function be used as a efficient change of variables?
    # The value of the function will probably have to be computed numerically
    # through a quadrature rule, but that should not be a problem. Test this?
end    

function derivative(f::Window{A,B,S},x) where {A,B,S}
    if  A ≤ x ≤ B
        return 1
    elseif 0 ≤ x < A
        u = (A - x) / (A)
        return exp(S*exp(-1/u)/(u-1))
    elseif B < x ≤ 1
        u = (x - B) / (1 - B)
        return exp(S*exp(-1/u)/(u-1))
    else
        return 0
    end
end    

"""
    struct Duffy{N} <: ChangeOfVariables
    
Change of variables mapping the reference `HyperCube{N}` to the `Simplex{N}`
with the property that the jacobian vanishes at the `𝐞₁` vertex of the simplex.
"""
struct Duffy{N} <: SingularityHandler{ReferenceTriangle} end

domain(::Duffy{2}) = ReferenceSquare()
range(::Duffy{2})  = ReferenceTriangle()

function (::Duffy{2})(u)
    SVector(u[1],(1-u[1])*u[2])
end    

function jacobian(::Duffy{2},u)
    SMatrix{2,2,Float64}(1,0,-u[2],(1-u[1]))
end    

struct TensorProductQuadratureHandler{S} <: SingularityHandler{ReferenceSquare} 
    cov::S
end

domain(::TensorProductQuadratureHandler) = ReferenceSquare()

function TensorProductQuadratureHandler(q...)
    TensorProductQuadratureHandler(q)
end    

function (f::TensorProductQuadratureHandler)(x)
    cov = f.cov    
    @assert length(cov) == length(x)
    svector(i->cov[i](x[i]),2)
end    

function jacobian(f::TensorProductQuadratureHandler,x)
    cov = f.cov    
    @assert length(cov) == length(x)
    jx = derivative(cov[1],x[1])
    jy = derivative(cov[2],x[2])
    SMatrix{2,2}(jx,0,0,jy)
end    
