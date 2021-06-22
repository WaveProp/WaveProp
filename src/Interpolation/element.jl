"""
    abstract type AbstractElement{D,T}

Fixed interpolation schemes over the domain `D<:AbstractReferenceShape`. This
means that the basis functions used for the interpolation are knonw from the
type. The type parameter `T` determines the [`return_type`](@ref) of the interpolant.

Instances `el` of `AbstractElement` are expected to implement:
- `el(x̂)`: evaluate the interpolation scheme at the coordinate `x̂ ∈ D`.
- `jacobian(el,x̂)` : evaluate the jacobian matrix of the interpolation at the
    coordinate `x ∈ D`.

!!! note
    For performance reasons, both `el(x̂)` and `jacobian(el,x̂)` should
    take as input a `StaticVector` and output a static vector or array.
"""
abstract type AbstractElement{D,T} end

function normal(el::AbstractElement, u)
    @assert u ∈ domain(el)
    jac = jacobian(el, u)
    normal(jac)
end

domain(::SType{<:AbstractElement{D}}) where {D <: AbstractReferenceShape} = D()

return_type(el::AbstractElement{D,T}) where {D,T} = T

domain_dimension(t::SType{<:AbstractElement}) = dimension(domain(t))

range_dimension(el::AbstractElement{R,T})        where {R,T} = length(T)

range_dimension(t::Type{<:AbstractElement{R,T}}) where {R,T} = length(T)

"""
    struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}

Standard element over `D <: AbstractReferenceShape` commonly used in finite
element methods. The underlying polynomial space is [`Pk{D,K}`](@ref), and it
evaluates to `vals::SVector{Np,T}` on the `reference_nodes` of the element.
"""
struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}
    vals::SVector{Np,T}
end

vals(el::LagrangeElement) = el.vals

# a contructor which infers the extra information from nodes
function LagrangeElement{R,K}(vals::SVector{Np,T}) where {R,K,T,Np}
    LagrangeElement{R,K,T}(vals)
end

# a contructor which infers the extra information from nodes
function LagrangeElement{D}(nodes::SVector{Np,T}) where {D,Np,T}
    LagrangeElement{D,Np,T}(nodes)
end

# constructor which converts each entry to a Point, and then creates an SVector
# of that.
function LagrangeElement{D}(vals) where {D}
    vals = SVector.(vals)
    LagrangeElement{D}(SVector(vals))
end

# a convenience constructor to allow things like LagrangeLine(a,b) instead of LagrangeLine((a,b))
function LagrangeElement{D}(nodes...) where {D}
    LagrangeElement{D}(nodes)
end

"""
    degree(el::LagrangeElement)

The polynomial degree of the element. A `LagrangeElement` of degree `K` and
domain `D` belongs to the space [`Pk{D,K}`](@ref).
"""
function degree(::SType{LagrangeElement{D,Np}})::Int where {D,Np}
    if D == ReferenceLine
        return Np-1
    elseif D == ReferenceTriangle
        K   = (-3 + sqrt(1+8*Np))/2
        return K
    elseif D == ReferenceTetrahedron
        notimplemented()
    elseif D == ReferenceSquare
        return sqrt(Np) - 1
    elseif D == ReferenceCube
        return Np^(1/3) - 1
    else
        notimplemented()
    end
end

function polynomial_space(::SType{LagrangeElement{D,Np}}) where {D,Np}
    K = degree(LagrangeElement{D,Np})
    Pk{D,K}()
end

"""
    const LagrangePoint = LagrangeElement{ReferencePoint}
"""
const LagrangePoint        = LagrangeElement{ReferencePoint}

"""
    const Point = LagrangePoint
"""
const Point = LagrangePoint

"""
    const LagrangeLine = LagrangeElement{ReferenceLine}
"""
const LagrangeLine        = LagrangeElement{ReferenceLine}

"""
    const LagrangeTriangle = LagrangeElement{ReferenceTriangle}
"""
const LagrangeTriangle    = LagrangeElement{ReferenceTriangle}

"""
    const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}
"""
const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}

"""
    const LagrangeSquare = LagrangeElement{ReferenceSquare}
"""
const LagrangeSquare   = LagrangeElement{ReferenceSquare}

"""
    reference_nodes(::LagrangeElement{D,Np,T})

Return the reference nodes on `D` used for the polynomial interpolation. The
function values on these nodes completely determine the interpolating
polynomial.

We use the same convention as `gmsh` for defining the reference nodes and their
order (see [node
ordering](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) on `gmsh`
documentation).
"""
function reference_nodes end

#=
Hardcode some basic elements.
TODO: Eventually this could/should be automated, at least for the LagrangeElements.
=#

# LagrangePoint represents the map () -> vals
function (el::LagrangePoint{1})()
    return vals(el)
end

# P1 for ReferenceLine
function (el::LagrangeLine{2})(u)
    @assert u ∈ ReferenceLine()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1]
end
function jacobian(el::LagrangeLine{2}, u)
    @assert u ∈ ReferenceLine()
    v = vals(el)
    return hcat(v[2] - v[1])
end

# P2 for ReferenceLine
function (el::LagrangeLine{3})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (4 * v[3] - 3 * v[1] - v[2]) * u[1]  + 2 * (v[2] + v[1] - 2 * v[3]) * u[1]^2
end
function jacobian(el::LagrangeLine{3}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(4 * v[3] - 3 * v[1] - v[2] + 4 * (v[2] + v[1] - 2 * v[3]) * u[1])
end

# P1 for ReferenceTriangle
function (el::LagrangeTriangle{3})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2]
    # v[1]*(1-u[1]-u[2]) + v[2]*(u[1]) + v[3]*u[2]
end
function jacobian(el::LagrangeTriangle{3}, u)
    @assert u ∈ domain(el)
    v   = vals(el)
    jac = hcat( v[2] - v[1],
                v[3] - v[1])
    return jac
end

# P2 for ReferenceTriangle
function (el::LagrangeElement{ReferenceTriangle,6})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    return (1+u[2]*(-3+2u[2])+u[1]*(-3+2u[1]+4u[2]))*v[1] +
           u[1]*(-v[2]+u[1]*(2v[2]-4v[4])+4v[4]+u[2]*(-4v[4]+4v[5]-4v[6])) +
           u[2]*(-v[3]+u[2]*(2v[3]-4v[6])+4v[6])
end
function jacobian(el::LagrangeElement{ReferenceTriangle,6}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        (-3+4u[1]+4u[2])*v[1] - v[2] + u[1]*(4v[2]-8v[4]) + 4v[4] + u[2]*(-4v[4]+4v[5]-4v[6]),
        (-3+4u[1]+4u[2])*v[1] - v[3] + u[2]*(4v[3]-8v[6]) + u[1]*(-4v[4]+4v[5]-4v[6]) + 4v[6]
    )
end

# P1 for ReferenceSquare
function (el::LagrangeElement{ReferenceSquare,4})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[4] - v[1]) * u[2] + (v[3] + v[1] - v[2] - v[4]) * u[1] * u[2]
end
function jacobian(el::LagrangeElement{ReferenceSquare,4}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        ((v[2] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[2]),
        ((v[4] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[1])
    )
end

# P1 for ReferenceTetrahedron
function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2] + (v[4] - v[1]) * u[3]
end
function jacobian(el::LagrangeElement{ReferenceTetrahedron,4}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        (v[2] - v[1]),
        (v[3] - v[1]),
        (v[4] - v[1])
    )
end
