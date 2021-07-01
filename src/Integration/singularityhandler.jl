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
image(::IMT)  = ReferenceLine()

function (f::IMT{A,P})(x) where {A,P}
    exp(A * (1 - 1 / x[1]^P))
end

derivative(f::IMT{A,P},x) where {A,P} = f(x) * A * P * 1 / x[1]^(P + 1)

jacobian(f::IMT,x) = derivative(f, x) |> SMatrix{1,1}

"""
    struct Kress{P} <: AbstractSingularityHandler{ReferenceLine}

Change of variables mapping `[0,1]` to `[0,1]` with the property that the first
`P` derivatives of the transformation vanish at `x=0`.
"""
struct Kress{P} <: AbstractSingularityHandler{ReferenceLine}
end
Kress(;order=5) = Kress{order}()

domain(k::Kress) = ReferenceLine()
image(k::Kress)  = ReferenceLine()

# NOTE: fastmath is needed here to allow for various algebraic simplifications
# which are not exact in floating arithmetic. Maybe reorder the operations *by
# hand* to avoid having to use fastmath? In any case, benchmark first.
function (f::Kress{P})(y) where {P}
    x = y[1]
    v = (x) -> (1 / P - 1 / 2) * ((1 - x))^3 + 1 / P * ((x - 1)) + 1 / 2
    return 2v(x)^P / (v(x)^P + v(2 - x)^P)
end

function derivative(f::Kress{P}, y) where {P}
    x = y[1]
    v = (x) -> (1 / P - 1 / 2) * ((1 - x))^3 + 1 / P * ((x - 1)) + 1 / 2
    vp = (x) -> -3 * (1 / P - 1 / 2) * ((1 - x))^2 + 1 / P
    return 2 * (P * v(x)^(P - 1) * vp(x) * (v(x)^P + v(2 - x)^P) - (P * v(x)^(P - 1) * vp(x) - P * v(2 - x)^(P - 1) * vp(2 - x) ) * v(x)^P ) /
        (v(x)^P + v(2 - x)^P)^2
end

jacobian(f::Kress,x) = derivative(f, x) |> SMatrix{1,1}

"""
    struct KressP{P} <: AbstractSingularityHandler{ReferenceLine}

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
image(k::KressP)  = ReferenceLine()

@fastmath function (f::KressP{P})(y) where {P}
    x = y[1]
    v = (x) -> (1 / P - 1 / 2) * ((1 - 2x))^3 + 1 / P * ((2x - 1)) + 1 / 2
    return v(x)^P / (v(x)^P + v(1 - x)^P)
end

@fastmath function derivative(f::KressP{P}, y) where {P}
    x = y[1]
    v =  (x) -> (1 / P - 1 / 2) * ((1 - 2x))^3 + 1 / P * ((2x - 1)) + 1 / 2
    vp = (x) -> -6 * (1 / P - 1 / 2) * ((1 - 2x))^2 + 2 / P
    return (P * v(x)^(P - 1) * vp(x) * (v(x)^P + v(1 - x)^P) - (P * v(x)^(P - 1) * vp(x) - P * v(1 - x)^(P - 1) * vp(1 - x) ) * v(x)^P ) /
        (v(x)^P + v(1 - x)^P)^2
end

jacobian(f::KressP,x) = derivative(f, x) |> SMatrix{1,1}

"""
    struct Duffy <: AbstractSingularityHandler{RefereceTriangle}

Change of variables mapping the `ReferenceSquare` to the `RefereceTriangle` with
the property that the jacobian vanishes at the `(1,0)` vertex of the triangle.

Useful for integrating functions with a singularity on the `(1,0)` edge of the
reference triangle.
"""
struct Duffy <: AbstractSingularityHandler{ReferenceTriangle} end

domain(::Duffy) = ReferenceSquare()
image(::Duffy)  = ReferenceTriangle()

function (::Duffy)(u)
    SVector(u[1], (1 - u[1]) * u[2])
end

function jacobian(::Duffy, u)
    SMatrix{2,2,Float64}(1, 0, -u[2][1], (1 - u[1][1]))
end

# TODO: generalize to `N` dimensions
"""
    struct TensorProductSingularityHandler{S} <: AbstractSingularityHandler{ReferenceSquare}

A tensor product of two one-dimensional `AbstractSingularityHandler`s for performing
integration over the `ReferenceSquare`.
"""
struct TensorProductSingularityHandler{S} <: AbstractSingularityHandler{ReferenceSquare}
    shandler::S
end

domain(::TensorProductSingularityHandler) = ReferenceSquare()
image(::TensorProductSingularityHandler) = ReferenceSquare()

function TensorProductSingularityHandler(q...)
    TensorProductSingularityHandler(q)
end

TensorProductSingularityHandler{Tuple{P,Q}}() where {P,Q} = TensorProductSingularityHandler(P(), Q())

function (f::TensorProductSingularityHandler)(x)
    shandler = f.shandler
    @assert length(shandler) == length(x)
    N = length(shandler)
    svector(i -> shandler[i](x[i]), N)
end

function jacobian(f::TensorProductSingularityHandler, x)
    shandler = f.shandler
    @assert length(shandler) == length(x)
    N = length(shandler)
    if N == 2
        jx = derivative(shandler[1], x[1])
        jy = derivative(shandler[2], x[2])
        return SMatrix{2,2}(jx, 0, 0, jy)
    else
        notimplemented()
    end
end
