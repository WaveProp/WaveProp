"""
    abstract type AbstractQuadratureRule{D<:AbstractReferenceShape}

A quadrature rule for integrating a function over the domain `D`.

An instance `q` of `AbstractQuadratureRule{D}` is expected to implement the
following methods:

- `q()` : return the nodes `x` and weights `w` of the quadrature rule on the
  reference domain `D`. For performance reasons, the result should depend only
  on the type of `q`.
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

## Define some one-dimensional quadrature rules

"""
    struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine}

`N`-point trapezoidal rule for integrating a function over the interval `[0,1]`.

For periodic functions over `[0,1]`, see [`TrapezoidalP`](@ref).

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
    h = 2 / (N - 1)
    x = [-1.0 + k * h for k in 0:N - 1]
    w = [h for _ in 1:N]
    w[1]   = h / 2
    w[end] = h / 2
    # convert to static arrays
    xs = svector(i -> SVector(0.5 * (x[i] + 1)), N)
    ws = svector(i -> w[i] / 2, N)
    return xs, ws
end

"""
    struct TrapezoidalP{N} <: AbstractQuadratureRule{ReferenceLine}

Open trapezoidal rule. Useful for periodic functions since it does not duplicate
the boundary nodes, void duplication.
"""
struct TrapezoidalP{N} <: AbstractQuadratureRule{ReferenceLine} end

TrapezoidalP(n::Int) = TrapezoidalP{n}()

# open trapezoidal rule for periodic functions
function _trapezoidalP(n)
    h = 1 / n
    x = [(k-0.5) * h for k in 1:n]
    w = [h for _ in 1:n]
    return x, w
end

function (q::TrapezoidalP{N})() where {N}
    x, w = _trapezoidalP(N)
    # convert to static arrays
    xs = SVector{N}(SVector{1}.(x))
    ws = SVector{N}(w)
    return xs, ws
end

"""
    struct Fejer{N}

`N`-point Fejer's first quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ N-1`.
"""
struct Fejer{N} <: AbstractQuadratureRule{ReferenceLine} end

Fejer(n::Int) = Fejer{n}()

# N point fejer quadrature integrates all polynomials up to degree N-1
order(::Fejer{N}) where {N} = N-1
function (q::Fejer{N})() where {N}
    theta = [(2j - 1) * π / (2 * N) for j = 1:N]
    x = -cos.(theta)
    w = zero(x)
    for j in 1:N
        tmp = 0.0
        for l in 1:floor(N / 2)
            tmp += 1 / (4 * l^2 - 1) * cos(2 * l * theta[j])
        end
        w[j] = 2 / N * (1 - 2 * tmp)
    end
    xs = svector(i -> SVector(0.5 * (x[i] + 1)), N)
    ws = svector(i -> w[i] / 2, N)
    return xs, ws
end

"""
    struct GaussLegendre{N}

`N`-point Gauss-Legendre quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ 2N-1`.
"""
struct GaussLegendre{N} <: AbstractQuadratureRule{ReferenceLine} end

"""
    GaussLegendre(;order)

Construct a `GaussLegendre` of the desired order over the `[0,1]` interval.
"""
function GaussLegendre(;order=p)
    N  = ceil(Int,(order + 1) /  2)
    GaussLegendre{N}()
end

GaussLegendre(n::Int) = GaussLegendre{n}()

# N point Gauss quadrature integrates all polynomials up to degree 2N-1, yielding
# an error of order 2N
order(q::GaussLegendre{N}) where {N} = 2*N-1

@generated function (q::GaussLegendre{N})() where {N}
    x, w  = gauss(N) # This is a quadgk function. Gives integral in [-1,1]. Converted to [0,1] below
    xs   = svector(i -> SVector(0.5 * (x[i] + 1)), N)
    ws   = svector(i -> w[i] / 2, N)
    return :($xs, $ws)
end

"""
    refine(q::AbstractQuadratureRule,[k=2])

Generate a similar quadrature rule, but with `k`-times as many quadrature nodes.
"""
function refine(q::GaussLegendre{N},k=2) where {N}
    GaussLegendre(Int(N*k))
end


"""
    struct Gauss{D,N} <: AbstractQuadratureRule{D}

Tabulated `N`-point symmetric Gauss quadrature rule for integration over `D`.
"""
struct Gauss{D,N} <: AbstractQuadratureRule{D}
    # gauss quadrature should be constructed using the order, and not the number
    # of nodes. This ensures you don't instantiate quadratures which are not
    # tabulated.
    function Gauss(;domain,order)
        if domain isa ReferenceTriangle
            msg = "quadrature of order $order not available for ReferenceTriangle"
            haskey(TRIANGLE_GAUSS_ORDER_TO_NPTS,order) || error(msg)
            n = TRIANGLE_GAUSS_ORDER_TO_NPTS[order]
        elseif domain isa ReferenceTetrahedron
            msg = "quadrature of order $order not available for ReferenceTetrahedron"
            haskey(TETRAHEDRON_GAUSS_ORDER_TO_NPTS,order) || error(msg)
            n = TETRAHEDRON_GAUSS_ORDER_TO_NPTS[order]
        else
            error("Tabulated Gauss quadratures only available for `ReferenceTriangle` or `ReferenceTetrahedron`")
        end
        return new{typeof(domain),n}()
    end
end

function order(q::Gauss{ReferenceTriangle,N}) where {N}
    TRIANGLE_GAUSS_NPTS_TO_ORDER[N]
end

function order(q::Gauss{ReferenceTetrahedron,N}) where {N}
    TETRAHEDRON_GAUSS_NPTS_TO_ORDER[N]
end

@generated function (q::Gauss{D,N})() where {D,N}
    x, w = _get_gauss_qnodes_and_qweights(D, N)
    return :($x, $w)
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
function integration_measure(jac::SMatrix)
    M,N = size(jac)
    if M == N
        det(jac) |> abs
    else
        transpose(jac)*jac |> det |> sqrt
    end
end

"""
    qrule_for_reference_shape(ref,order)

Given a `ref`erence shape and a desired quadrature `order`, return
an appropiate quadrature rule.
"""
function qrule_for_reference_shape(ref,order)
    if ref isa ReferenceLine
        return GaussLegendre(;order)
    elseif ref isa ReferenceSquare
        qx = qrule_for_reference_shape(ReferenceLine(),order)
        qy = qx
        return TensorProductQuadrature(qx,qy)
    elseif ref isa ReferenceTriangle
        return Gauss(;domain=ref,order=order)
    elseif ref isa ReferenceTetrahedron
        return Gauss(;domain=ref,order=order)
    else
        error("no appropriate quadrature rule found.")
    end
end
