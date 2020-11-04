"""
    abstract type AbstractQuadratureRule{D<:AbstractReferenceShape}
    
A quadrature rule for integrating a function over the domain `D`.

An instance `q` of `AbstractQuadratureRule{D}` is expected to implement the
following methods:

- `q()` : return the nodes `x` and weights `w` of the quadrature rule on the
  reference domain `D`. For performance reasons, the result shoudl depend only
  on the type of `q`. 
- `q(el)` : return the nodes `x` and weights `w` of the quadrature rule on the
  elemenent `D`. This assumes that `domain(q)==domain(el)`, so that the element
  quadrature can be computed by *pushing forward* a reference quadrature to `el`.  
"""
abstract type AbstractQuadratureRule{D} end

domain(q::AbstractQuadratureRule{D}) where {D} = D()

"""
    (q::AbstractQuadratureRule)()

Return the quadrature nodes `x` and weights `w` on the `domain(q)`.
"""
function (q::AbstractQuadratureRule) end

"""
    (q::AbstractQuadratureRule)(el)

Return the quadrature nodes `x` and weights `w` for integrating over `el`. Here
`el` can represent an element, or a change of variables, as long as
`domain(el)==domain(q)`. 

The *lifted* quadrature is computed by mapping the reference quadrature through
`el`. This requires `el` to support the methods `el(x̂)` and `jacobian(el,x̂)`.
"""
function (q::AbstractQuadratureRule)(el)
    @assert domain(el) == domain(q) "the domains of the `q` and `el` must agree"
    x̂,ŵ = q()
    x   = map(x->el(x),x̂)
    w   = map(zip(x̂,ŵ)) do (x̂,ŵ)
        μ = measure(el,x̂)
        μ*prod(ŵ)
    end 
    return x,w
end

"""
    integrate(f,q::AbstractQuadrature)
    integrate(f,x,w)

Integrate the function `f` using the quadrature rule `q`. This is simply
`sum(f.(x) .* w)`, where `x` and `w` are the quadrature nodes and weights, respectively.
"""
function integrate(f,q::AbstractQuadratureRule)
    x,w = q()
    integrate(f,x,w)
end    

function integrate(f,x,w)
    sum(zip(x,w)) do (x,w)
        f(x)*prod(w)
    end
end    

# overload quadgk for integration over reference shapes. Useful for various
# testing purposes.
"""
    integrate(f,s::AbstractReferenceShape)

Use `quadgk` to (adaptively) integrate a function over the reference shape `s`. 

This is used mostly for testing purposes. It returns only the value of the
integral (i.e. without the error estimate).
"""
integrate(f,l::AbstractReferenceShape) = quadgk(f,typeof(l))

integrate(f,::Type{ReferenceLine})     = quadgk(f,0,1)[1]

function integrate(f,::Type{ReferenceSquare})
    I    = x-> quadgk(y->f(Point(x,y)),0,1)[1]
    out  = quadgk(I,0,1)[1]
end    

function integrate(f,::Type{ReferenceTriangle})
    I    = x -> quadgk(y->f(Point(x,y)),0,1-x)[1]
    out  = quadgk(I,0,1)[1]
end

## Define some one-dimensional quadrature rules

"""
    struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine}

`N`-point trapezoidal rule for integrating a function over the interval `[0,1]`.

Note that for analytic 1-periodic functions, this rule will converge
exponentially fast. 
    
# Examples:
```julia
q    = Trapezoidal(10)
f(x) = exp(x)*cos(x)
integrate(f,q)
```
"""
struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine} end

Trapezoidal(n::Int) = Trapezoidal{n}()

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
    struct Fejer{N}

`N`-point Fejer's first quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ N-1`.
"""
struct Fejer{N} <: AbstractQuadratureRule{ReferenceLine} end

Fejer(n::Int) = Fejer{n}()

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
    struct GaussLegendre{N}

`N`-point Gauss-Legendre quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ 2N-1`.
"""
struct GaussLegendre{N} <: AbstractQuadratureRule{ReferenceLine} end

GaussLegendre(n::Int) = GaussLegendre{n}()

function (q::GaussLegendre{N})() where {N}
    x,w  = gauss(N) # This is a quadgk function. Gives integral in [-1,1]. Converted to [0,1] below
    xs   = svector(i->(0.5*(x[i]+1)),N) 
    ws   = svector(i->w[i]/2,N)
    return xs,ws
end
 
"""
    struct Gauss{D,N} <: AbstractQuadratureRule{D}
    
Tabulated `N`-point symmetric Gauss quadrature rule for integration over `D`.

This is currently implemented *by hand* for low values of `N` on triangles and tetrahedrons.
"""
struct Gauss{D,N} <: AbstractQuadratureRule{D} end

Gauss(ref,n) = Gauss{typeof(ref),n}()
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

"""
    TensorProductQuadrature{Q}

A tensor-product of one-dimension quadrature rules. Integrates over `[0,1]^d`,
where `d=length(quad)`.

# Examples
```julia
qx = Fejer(10)
qy = GaussLegendre(15)
q  = TensorProductQuadrature(qx,qy)
```
"""
struct TensorProductQuadrature{Q} <: AbstractQuadratureRule{ReferenceSquare} 
    quad::Q
end

function TensorProductQuadrature(q...)
    TensorProductQuadrature(q)
end    

# FIXME: instead of returning an iterator, the tensor product rule is currently
# returning the actual matrices. 
function (q::TensorProductQuadrature)()
    nodes   = map(q->q()[1],q.quad)    
    weights = map(q->q()[2],q.quad)    
    x = Iterators.product(nodes...) 
    w = Iterators.product(weights...) 
    return Point.(x), prod.(collect(w))
end    
