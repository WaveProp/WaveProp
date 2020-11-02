"""
    abstract type AbstractQuadratureRule{D<:AbstractReferenceShape}
    
A quadrature rule for integrating a function over the region `D`.
"""
abstract type AbstractQuadratureRule{D} end

domain(q::AbstractQuadratureRule{D}) where {D} = D()

function integrate(f,q::AbstractQuadratureRule)
    x,w = q()
    integrate(f,x,w)
end    

function integrate(f,x,w)
    mapreduce(+,zip(x,w)) do (x,w)
        f(x)*prod(w)
    end
end    

# One-dimensional quadratures

"""
    struct Trapezoidal{N}()

`N`-point trapezoidal rule for integrating a function over the interval `[0,1]`.

Note that for analytic periodic functions, the trapezoidal rule converges exponentially fast. 
"""
struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine} end

function (q::Trapezoidal{N})() where {N}
    h = 2/(N-1)
    x = [-1.0 + k*h for k in 0:N-1]
    w = [h for k in 1:N]
    w[1]   = h/2
    w[end] = h/2
    # convert to static arrays
    xs = svector(i->(0.5*(x[i]+1)),N) 
    ws = svector(i->w[i]/2,N)
    return xs,ws
end    

"""
    struct Fejer{N}()

`N`-point Fejer's first quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials up to degree `N-1`.
"""
struct Fejer{N} <: AbstractQuadratureRule{ReferenceLine} end

function (q::Fejer{N})() where {N}
    theta = [(2j-1)*π/(2*N) for j=1:N]
    x = -cos.(theta)
    w = zero(x)
    for j in 1:N
        tmp = 0.0
        for l in 1:floor(N/2)
            tmp += 1/(4*l^2-1) *cos(2*l*theta[j])
        end
        w[j] = 2/N * (1 - 2*tmp)
    end
    xs = svector(i->(0.5*(x[i]+1)),N) 
    ws = svector(i->w[i]/2,N)
    return xs,ws
end

"""
    struct GaussLegendre{N}()

`N`-point Gauss-Legendre quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials up to degree `2N-1`.
"""
struct GaussLegendre{N} <: AbstractQuadratureRule{ReferenceLine} end

function (q::GaussLegendre{N})() where {N}
    x,w  = gauss(N) # gives integral in [-1,1]. Converted to [0,1] below
    xs   = svector(i->(0.5*(x[i]+1)),N) 
    ws   = svector(i->w[i]/2,N)
    return xs,ws
end

# struct DoubleExponential{N} <: AbstractQuadratureRule{N} end

# """
#     double_exponential(n::Integer) -> (x,w)

# Compute nodes `x` and weights `w` of the double expoenential quadrature rule for integrating a function over `[-1,1]`.
# """
# function double_exponential(n::Integer,h=_estimate_h(n))
#     @assert isodd(n) "double exponential requires odd number of points."
#     x̂ = (t) -> tanh(π/2*sinh(t))
#     ŵ = (t) -> h*π/2*cosh(t)/cosh(π/2*sinh(t))^2
#     krange = -(n-1)÷2 : (n-1)÷2
#     x = [x̂(k*h) for k in krange]
#     w = [ŵ(k*h) for k in krange]
#     return x,w
# end

# # FIXME: find out a how to estimate the stepsize h in the double exponential
# # formula. I think what you want is to find `h` s.t. n*h = constant, where the
# # constant should depend on the desired precision. 
# function _estimate_h(n;tol=1e-16)
#     x̂ = (t) -> tanh(π/2*sinh(t))
#     ŵ = (t) -> h*π/2*cosh(t)/cosh(π/2*sinh(t))^2
#     h = 2^(-8)
#     while true
#         if ŵ(n*h) < tol || 1-x̂(n*h) < tol
#             return h    
#         else
#             h = 2*h
#         end    
#     end
# end    

"""
    struct Gauss{D,N} <: AbstractQuadratureRule{D}
    
Tabulated `N`-point symmetric Gauss quadrature rule for integration over `D`.

This is currently implemented for low values of `N` on triangles and tetrahedrons.
"""
struct Gauss{D,N} <: AbstractQuadratureRule{D} end

Gauss(ref;n) = Gauss{typeof(ref),n}()

function (q::Gauss{ReferenceTriangle,N})() where {N}
    if N == 1
        x = svector(i->Point(1/3,1/3),1)
        w = svector(i->1/2,1)
    elseif N == 3
        x = SVector(Point(1/6,1/6),
                    Point(2/3,1/6),
                    Point(1/6,2/3))
        w = svector(i->1/6,3)
    else
        @notimplemented
    end            
    return x,w
end

function get_order(q::Gauss{ReferenceTriangle,N}) where {N}
    if N == 1
        return 1
    elseif N == 3
        return 3
    else    
        @notimplemented        
    end    
end    

function (q::Gauss{ReferenceTetrahedron,N})() where {N}
    if N == 1
        x = SVector((Point(1/4,1/4,1/4),))
        w = SVector(1/6)
    elseif N == 4
        a = (5-√5)/20
        b = (5+3*√5)/20
        x = SVector(Point(a,a,a),
                    Point(a,a,b),
                    Point(a,b,a),
                    Point(b,a,a)
                )
        w = svector(i->1/24,4)
    else
        @notimplemented
    end            
    return x,w
end

struct TensorProduct{Q} <: AbstractQuadratureRule{ReferenceSquare} 
    quad::Q
end

function TensorProduct(q...)
    TensorProduct(q)
end    

# FIXME: instead of returning an iterator, the tensor product rule is currently
# returning the actual matrices. 
function (q::TensorProduct)()
    nodes   = map(q->q()[1],q.quad)    
    weights = map(q->q()[2],q.quad)    
    x = Iterators.product(nodes...) 
    w = Iterators.product(weights...) 
    return Point.(x), prod.(collect(w))
end    

function (qrule::AbstractQuadratureRule)(el)
  x̂,ŵ = qrule()
  x,w = _push_forward_quad(el,x̂,ŵ)
  return x,w
end

function _push_forward_quad(cov,x̂,ŵ)
    x   = map(x->cov(x),x̂)
    w   = map(zip(x̂,ŵ)) do (x̂,ŵ)
        μ = measure(cov,x̂)
        μ*prod(ŵ)
    end 
    return x,w
end

function push_forward_quad_with_normal(el,qrule)
    @assert domain(el) == domain(qrule)    
    x̂,ŵ = qrule()
    _push_forward_quad_with_normal(el,x̂,ŵ)
end

function _push_forward_quad_with_normal(el,x̂,ŵ)
    x,w = _push_forward_quad(el,x̂,ŵ)
    ν   = map(x->normal(el,x),x̂)
    return x,w,ν
end
