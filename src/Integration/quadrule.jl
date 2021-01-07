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
    qnodes(Y)

Return the quadrature nodes associated with `Y`.
"""
qnodes(q::AbstractQuadratureRule) = q()[1]

"""
    qweights(Y)

Return the quadrature weights associated with `Y`.
"""
qweights(q::AbstractQuadratureRule) = q()[2]

"""
    qnormals(Y)

Return the normal vector at the quadrature nodes of `Y.
"""
qnormals(q::AbstractQuadratureRule) = error("abstract quadrature rule has no normal")

"""
    (q::AbstractQuadratureRule)()

Return the quadrature nodes `x` and weights `w` on the `domain(q)`.
"""
function (q::AbstractQuadratureRule)() end

"""
    integrate(f,q::AbstractQuadrature)
    integrate(f,x,w)

Integrate the function `f` using the quadrature rule `q`. This is simply
`sum(f.(x) .* w)`, where `x` and `w` are the quadrature nodes and weights, respectively.
"""
function integrate(f,q::AbstractQuadratureRule)
    x,w = q()
    if domain(q) == ReferenceLine()
        return integrate(x->f(x[1]),x,w)
    else
        return integrate(f,x,w)
    end    
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
integrate(f,l::AbstractReferenceShape;kwargs...)     = integrate(f,typeof(l);kwargs...)

integrate(f,::Type{ReferenceLine};kwargs...) = quadgk(f,0,1;kwargs...)[1]

function integrate(f,::Type{ReferenceSquare})
    I    = x-> quadgk(y->f(SVector(x,y)),0,1)[1]
    out  = quadgk(I,0,1)[1]
end    

# hacky way to compute integration over reference triangles. Only used for
# testing purposes to avoid having to compute integrals analyically when testing.
function integrate(f,::Type{ReferenceTriangle})
    I    = x -> quadgk(y->f(SVector(x,y)),0,1-x)[1]
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
    xs = svector(i->SVector(0.5*(x[i]+1)),N) 
    ws = svector(i->w[i]/2,N)
    return xs,ws
end   

"""
    struct TrapezoidalP{N} <: AbstractQuadratureRule{ReferenceLine}

Like [`Trapezoidal{N}`](@ref), but assumes the function is 1-periodic, and so the
nodes and weights at the endSVectors of the interval `[0,1]` need not be duplicated.
"""
struct TrapezoidalP{N} <: AbstractQuadratureRule{ReferenceLine} end

TrapezoidalP(n::Int) = TrapezoidalP{n}()

# open periodic trapezoidal rule on [0,1]
function _trapezoidalP(n)    
    h = 1/n
    x = [(k-0.5)*h for k in 1:n]
    w = [h   for k in 1:n]
    return x,w
end    

function (q::TrapezoidalP{N})() where {N}
    x,w = _trapezoidalP(N)
    # convert to static arrays
    xs = SVector{N}(SVector{1}.(x))
    ws = SVector{N}(w)
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
    xs = svector(i->SVector(0.5*(x[i]+1)),N) 
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

@generated function (q::GaussLegendre{N})() where {N}
    x,w  = gauss(N) # This is a quadgk function. Gives integral in [-1,1]. Converted to [0,1] below
    xs   = svector(i->SVector(0.5*(x[i]+1)),N) 
    ws   = svector(i->w[i]/2,N)
    return :($xs,$ws)
end

function refine(q::GaussLegendre{N},k=2) where {N}
    GaussLegendre(Int(N*k))
end    

"""
    struct Gauss{D,N} <: AbstractQuadratureRule{D}
    
Tabulated `N`-point symmetric Gauss quadrature rule for integration over `D`.

This is currently implemented *by hand* for low values of `N` on triangles and tetrahedrons.
"""
struct Gauss{D,N} <: AbstractQuadratureRule{D} end

Gauss(ref,n) = Gauss{typeof(ref),n}()
Gauss(ref;n) = Gauss{typeof(ref),n}()

@generated function (q::Gauss{ReferenceTriangle,N})() where {N}
    if N == 1
        x = svector(i->SVector(1/3,1/3),1)
        w = svector(i->1/2,1)
    elseif N == 3
        x = SVector(SVector(1/6,1/6),
                    SVector(2/3,1/6),
                    SVector(1/6,2/3))
        w = svector(i->1/6,3)
    else
        notimplemented()
    end            
    return :($x,$w)
end

@generated function (q::Gauss{ReferenceTetrahedron,N})() where {N}
    if N == 1
        x = SVector((SVector(1/4,1/4,1/4),))
        w = SVector(1/6)
    elseif N == 4
        a = (5-√5)/20
        b = (5+3*√5)/20
        x = SVector(SVector(a,a,a),
                    SVector(a,a,b),
                    SVector(a,b,a),
                    SVector(b,a,a)
                )
        w = svector(i->1/24,4)
    else
        notimplemented()
    end            
    return :($x,$w)
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

# FIXME: this is a workaround the need to easily construct a tensor quadrature
# based only on the the types of the quadratures. Useful in generated functions,
# but there is probably a better way
function TensorProductQuadrature{Tuple{Q1,Q2}}() where {Q1,Q2}
    TensorProductQuadrature(Q1(),Q2())
end    

function TensorProductQuadrature(q...)
    TensorProductQuadrature(q)
end    

# FIXME: the current implementation is rather obscure. How should we handle the
# product quadrature rules in general? Also make this into a generated function.
function (q::TensorProductQuadrature)()
    N       = length(q.quad)    
    nodes   = map(q->q()[1],q.quad)    
    weights = map(q->q()[2],q.quad)  
    nodes_iter   = Iterators.product(nodes...)
    weights_iter = Iterators.product(weights...)
    x = map(x->vcat(x...),nodes_iter)
    w = map(w->prod(w),weights_iter)
    return SArray(x), w
end    

# some helper functions

"""
    _qrule_for_reference_shape(ref,order)

Given a `ref`erence shape and a desired quadrature `order`, return 
an appropiate quadrature rule.
"""
function _qrule_for_reference_shape(ref,order)
    if ref isa ReferenceLine
        n = ((order + 1) ÷  2) + 1
        return GaussLegendre{n}()
    elseif ref isa ReferenceSquare
        n  = (order + 1)/2 |> ceil
        qx = GaussLegendre{n}()
        qy = qx
        return TensorProductQuadrature(qx,qy)
    elseif ref isa ReferenceTriangle
        if order <= 1
            return Gauss(ref,n=1) 
        elseif order <=2
            return Gauss(ref,n=3)     
        end
    elseif ref isa ReferenceTetrahedron
        if order <= 1
            return Gauss(ref;n=1) 
        elseif order <=2
            return Gauss(ref;n=4)
        end
    end    
    error("no appropriate quadrature rule found.")
end   

"""
    _qrule_for_element(E,p)

Given an element type `E`, return an appropriate quadrature of order `p`.
"""
function _qrule_for_element(E,order)
    _qrule_for_reference_shape(domain(E),order)
end    
