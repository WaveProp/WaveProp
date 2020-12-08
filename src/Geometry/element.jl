"""
    abstract type AbstractReferenceShape{N}
    
A reference polygon in `ℜᴺ`.

Reference shapes are used mostly for defining [`AbstractElement`](@ref)s as
transformations mapping an `AbstractReferenceShape` into some region of `ℜᴹ`. 

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for some
examples of concrete subtypes.
"""
abstract type AbstractReferenceShape{N} end

ambient_dimension(::AbstractReferenceShape{N}) where {N} = N
geometric_dimension(::AbstractReferenceShape{N}) where {N} = N

"""
    struct ReferenceLine
    
Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: AbstractReferenceShape{1}
end
Base.in(x,::ReferenceLine) = 0 ≤ x[1] ≤ 1
center(::Type{ReferenceLine}) = 0.5
center(::ReferenceLine)       = 0.5

"""
    struct ReferenceTriangle
    
Singleton type representing the triangle with vertices `(0,0),(0,1),(1,0)`
"""
struct ReferenceTriangle <: AbstractReferenceShape{2} 
end    
Base.in(x,::ReferenceTriangle) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1]

"""
    struct ReferenceSquare
    
Singleton type representing the square with vertices `(0,0),(0,1),(1,1),(1,0)`
"""
struct ReferenceSquare <: AbstractReferenceShape{2} 
end    
Base.in(x,::ReferenceSquare) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1
center(::Type{ReferenceSquare}) = SVector(0.5,0.5)
center(::ReferenceSquare)       = SVector(0.5,0.5)

"""
    struct ReferenceTetrahedron
    
Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: AbstractReferenceShape{3} 
end    
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1] && 0 ≤ x[3] ≤ 1 - x[1] - x[2]

# TODO: generalize structs above to `ReferenceSimplex{N}` and `ReferenceCuboid{N}`

"""
    abstract type AbstractElement{D,N}

Abstract shape given by the image of a parametrization with domain
`D<:AbstractReferenceShape`.

The type parameter `N` represents the dimension of the embedding space. The
geometric dimension of the element can be obtained from the geometric dimension
of its domain `D`. 

Instances `el` of `AbstractElement` are expected to implement.

- `el(x̂)`: evaluate the parametrization defining the element at the parametric
    coordinates `x̂ ∈ D`.
- `jacobian(el,x̂)` : evaluate the jacobian matrix of the parametrization at the
    parametric coordinate `x ∈ D`. For performance reasons, this should return
    an `SMatrix` of size `M × N`, where `M` is the [`ambient_dimension`](@ref)
    of `el` and `N` is the [`geometric_dimension`](@ref) of `el`, respectively.
"""
abstract type AbstractElement{D,N} end

const AbstractLine{N} = AbstractElement{ReferenceLine,N}
const AbstractTriangle{N} = AbstractElement{ReferenceTriangle,N}

"""
    (el::AbstractElement)(x)

Evaluate the underlying parametrization of the element `el` at point `x`. This is the push-forward map for the element. 
"""
function (el::AbstractElement)(x) 
    abstractmethod(typeof(el))
end

"""
    domain(el::AbstractElement)

Return an instance of the singleton type `R`; i.e. the reference element.
"""
domain(::Type{<:AbstractElement{R}}) where {R<:AbstractReferenceShape} = R()
domain(el::AbstractElement) = domain(typeof(el))

"""
    measure(τ,u)

The integration measure `μ` of the transformation `τ` so that 
```math
\\int_\\tau f(y) ds_y = \\int_{\\hat{\\tau}} f(\\tau(\\hat{y})) \\mu(\\hat{y}) d\\hat{y}
```
where `` \\hat{\\tau} `` is the reference element.
"""
function measure(el,u)
    jac = jacobian(el,u)
    g   = transpose(jac)*jac |> det
    μ   = sqrt(g)
    return μ
end    

"""
    normal(el::AbstractElement,x̂)

The outer normal vector for the `el` at the parametric coordinate `x̂ ∈
domain(el)`.
    
Note that `normal` is only defined for co-dimension one elements. 
"""
function normal(el::AbstractElement,u)
    N = ambient_dimension(el)
    M = geometric_dimension(el)
    msg = "computing the normal vector requires the element to be of co-dimension one."
    @assert N-M == 1 msg
    if M == 1 # a line in 2d
        t⃗ = jacobian(el,u)
        n⃗ = SVector(t⃗[2],-t⃗[1])
        return n⃗/norm(n⃗)
    elseif M == 2 # a surface in 3d
        j  = jacobian(el,u)    
        t⃗₁ = j[:,1]
        t⃗₂ = j[:,2]
        n⃗  = cross(t⃗₁,t⃗₂)
        return n⃗/norm(n⃗)
    else
        notimplemented()    
    end            
end    

"""
    geometric_dimension(el::AbstractElement{R})

Return the geometric dimension of `el`, i.e. the number of variables needed to locally
parametrize the element.
"""
geometric_dimension(t::Type{<:AbstractElement}) = geometric_dimension(domain(t))
geometric_dimension(el::AbstractElement)        = geometric_dimension(typeof(el))

"""
    ambient_dimension(el::AbstractElement)

Return the dimension of the ambient space where `el` lives.
"""
ambient_dimension(el::AbstractElement{R,N})   where {R,N} = N
ambient_dimension(t::Type{<:AbstractElement{R,N}}) where {R,N} = N

boundary(el::AbstractLine) = el(0),el(1)

function jacobian(el::AbstractLine,u;h=sqrt(eps()))
    (el(u+h)-el(u))/h
end    

"""
    abstract type PolynomialElement{R,M} <: AbstractElement{R,M}
    
Abstract elements whose parametrization is a polynomial.
"""
abstract type PolynomialElement{R,N} <: AbstractElement{R,N} end

"""
    struct LagrangeElement{R,Np,N,T}        
    
# Fields:
- `nodes::SVector{Np,Point{N,T}}`

A lagrange element is represented as a polynomial mapping the `Np` reference
lagrangian nodes of the reference element `R` into `nodes`.

The element's parametrization is fully determined by the image of the `Np`
reference points, which are stored in `nodes`.
"""
struct LagrangeElement{R,Np,N,T} <: PolynomialElement{R,N}
    nodes::SVector{Np,Point{N,T}}
end

# a contructor which infers the extra information from nodes
function LagrangeElement{R}(nodes::SVector{Np,Point{N,T}}) where {R,Np,N,T} 
    LagrangeElement{R,Np,N,T}(nodes)
end

# constructor which converts each entry to a Point, and then creates an SVector
# of that.
function LagrangeElement{R}(nodes...) where {R} 
    nodes = Point.(nodes)
    LagrangeElement{R}(SVector(nodes))
end

# define some aliases for convenience
const LagrangeLine        = LagrangeElement{ReferenceLine}
const LagrangeTriangle    = LagrangeElement{ReferenceTriangle}
const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}
const LagrangeRectangle   = LagrangeElement{ReferenceSquare}

"""
    line(a,b)

Create a straight line connecting points `a` and `b`. This returns an instance
of [`LagrangeLine`](@ref).
"""
line(a,b) = LagrangeLine(a,b)

"""
    triangle(a,b,c)

Create a triangle connecting points `a`, `b`, and `c`. This returns an instance
of [`LagrangeTriangle`](@ref).
"""
triangle(a,b,c) = LagrangeTriangle(a,b,c)

"""
    tetrahedron(a,b,c,d)

Create a tetrahedron with vertices `a`, `b`, `c` and `d`. This returns an instance
of [`LagrangeTetrahedron`](@ref).
"""
tetrahedron(a,b,c,d) = LagrangeTetrahedron(a,b,c,d)

"""
    rectangle(a,b,c,d)

Create a rectangle with vertices `a`, `b`, `c` and `d`. This returns an instance
of [`LagrangeRectangle`](@ref).
"""
rectangle(a,b,c,d) = LagrangeRectangle(a,b,c,d)

nodes(el::LagrangeElement) = el.nodes

get_order(el::LagrangeElement{ReferenceLine,Np}) where {Np} = Np + 1

# For a lagrangian element of order P on a triangle, there are Np = (P+1)(P+2)/2
# elements
function get_order(el::LagrangeElement{ReferenceTriangle,Np}) where {Np} 
    # TODO: make this generic by solving quadratic equation and assert integer solution
    if Np == 3
        return 1
    elseif Np == 6 
        return 2
    else    
        msh = "unable to determine order for element of type $(el) containing Np=$(Np)        interpolation points. Make sure Np = (P+1)(P+2)/2 for some integer P."
        return error(msh)        
    end    
end 

"""
    get_reference_nodes(::LagrangeElement{R,Np,N,T})

Return the nodes on the reference element for a given `el::LagrangeElement`. 

For performance reasons, this should depend solely on the type of `LagrangeElement` so that no runtime computations need to be performed. The method returns a `SVector` with `Np` elements of type `Point{M,T}` where `M` is the dimension of the reference shape `R`.
"""
function get_reference_nodes end

function get_reference_nodes(el::LagrangeElement{ReferenceLine,Np}) where {Np <: Integer}
    dx = 1 / (Np - 1)
    svector(i -> Point((i - 1) * dx), Np)
end

function get_reference_nodes(el::LagrangeElement{ReferenceTriangle,Np}) where {Np <: Integer}
    if Np == 3
        return SVector(Point(0., 0.), Point(1., 0.), Point(1.0, 0.0))    
    elseif Np == 6
        return SVector(
                        Point(0., 0.),Point(0.5, 0.0),Point(0., 1.0),
                        Point(0.0, 0.5),Point(0.5, 0.5),
                        Point(0.0, 1.0)
                )        
    else    
        notimplemented()
    end
end

function (el::LagrangeLine{2})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    x = nodes(el)
    x[1] + (x[2] - x[1]) * u[1]    
end    

function (el::LagrangeLine{3})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    x = nodes(el)
    x[1] + (4*x[3] - 3*x[1] - x[2])*u[1]  + 2*(x[2]+x[1]-2*x[3])*u[1]^2
end    

function jacobian(el::LagrangeLine{3},u)
    N = ambient_dimension(el)        
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    x = nodes(el)
    SMatrix{N,1}((4*x[3] - 3*x[1] - x[2] + 4*(x[2]+x[1]-2*x[3])*u[1])...)
end    

function (el::LagrangeElement{ReferenceTriangle,3})(u)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle()
    x = nodes(el)
    x[1] + (x[2] - x[1]) * u[1] + (x[3] - x[1]) * u[2]
end 

function (el::LagrangeElement{ReferenceSquare,4})(u)
    @assert length(u) == 2    
    @assert u ∈ ReferenceSquare()
    x = nodes(el)
    x[1] + (x[2] - x[1]) * u[1] + (x[4] - x[1]) * u[2] + (x[3]+x[1]-x[2]-x[4]) * u[1]*u[2]
end 

function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert length(u) == 3    
    @assert u ∈ ReferenceTetrahedron()
    x = nodes(el)
    x[1] + (x[2] - x[1]) * u[1] + (x[3] - x[1]) * u[2] + (x[4] - x[1]) * u[3]
end 

"""
    jacobian(el::LagrangeElement,x)

Evaluate the jacobian of the underlying parametrization of the element `el` at point `x`. 
"""
function jacobian(el::LagrangeElement{ReferenceLine,2}, u)
    N = ambient_dimension(el)    
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    x = nodes(el)
    return SMatrix{N,1}((x[2] - x[1])...)
end    

function jacobian(el::LagrangeElement{ReferenceTriangle,3}, u) 
    N = ambient_dimension(el)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle()
    x = nodes(el)
    SMatrix{N,2}(
        (x[2] - x[1])..., 
        (x[3] - x[1])...
    )
end 

function jacobian(el::LagrangeElement{ReferenceSquare,4},u)
    N = ambient_dimension(el)
    @assert length(u) == 2    
    @assert u ∈ ReferenceSquare()
    x = nodes(el)
    hcat(
        ((x[2] - x[1]) + (x[3]+x[1]-x[2]-x[4])*u[2]),
        ((x[4] - x[1]) + (x[3]+x[1]-x[2]-x[4])*u[1])
    )
end 

function jacobian(el::LagrangeElement{ReferenceTetrahedron,4}, u)
    N = ambient_dimension(el)    
    @assert length(u) == 3   
    @assert u ∈ ReferenceTriangle()
    x = nodes(el)
    SMatrix{N,3}( 
        (x[2] - x[1])...,
        (x[3] - x[1])...,
        (x[4] - x[1])...
    )
end 

"""
    abstract type ParametricElement{R,N} <: AbstractElement{R,N}
    
An element given by some explicit parametric function.        
"""
abstract type AbstractParametricElement{R,N} <: AbstractElement{R,N} end

"""
    ParametricElement{D,N,F}

An element represented as the explicit mapping `f::F` with domain `D`.  

Unlike e.g. `LagrangeElements`, whose domain are `<:AbstractReferenceShape`s, a
`ParametricElement` has a domain given by a `HyperRectangle`. This means 
"""
struct ParametricElement{D,T,F} <: AbstractParametricElement{D,T}
    parametrization::F
    domain::D    
end

Base.eltype(p::ParametricElement{D,T,F}) where {D,T,F} = T

# The domain of a parametric element is the reference domain which can be used
# to describe it.
function domain(p::Type{ParametricElement{D,T,F}}) where {D,T,F}
    if D <: HyperRectangle{1}
        return ReferenceLine()    
    elseif D <: HyperRectangle{2}
        return ReferenceSquare()    
    else
        notimplemented()
    end        
end    
domain(p::ParametricElement) = domain(typeof(p))

geometric_dimension(p::ParametricElement) = geometric_dimension(domain(p))
ambient_dimension(p::ParametricElement)   = length(eltype(p))

# constructor which infers the return type of f by applying it to a point in the
# reference domain.
function ParametricElement(f,d)
    F = typeof(f)
    D = typeof(d)
    T = f(center(d)) |> typeof
    return ParametricElement{D,T,F}(f,d)
end  

function (el::ParametricElement)(u) 
    @assert u ∈ domain(el)
    rec = el.domain
    lc  = low_corner(rec)
    hc  = high_corner(rec)
    N = geometric_dimension(el)
    # map from reference domain to domain of element (confusing...)
    v = svector(N) do dim
        lc[dim] + (hc[dim]-lc[dim])*u[dim]
    end    
    el.parametrization(v)
end

function jacobian(el::ParametricElement,u::SVector) 
    @assert u ∈ domain(el)
    rec = el.domain
    lc  = low_corner(rec)
    hc  = high_corner(rec)
    N = geometric_dimension(el)
    # map from reference domain to domain of element (confusing...)
    v = svector(N) do dim
        lc[dim] + (hc[dim]-lc[dim])*u[dim]
    end        
    scal = svector(N) do dim
        (hc[dim]-lc[dim])
    end        
    ForwardDiff.jacobian(el.parametrization,v) * SDiagonal(scal)
end
jacobian(psurf::ParametricElement,s) = jacobian(psurf,SVector(s))

# define some aliases for convenience
const ParametricLine  = ParametricElement{HyperRectangle{1,<:Number}}

derivative(l::ParametricLine,s)  = ForwardDiff.derivative(l,s)
derivative2(l::ParametricLine,s) = ForwardDiff.derivative(s -> derivative(l,s),s)
derivative3(l::ParametricLine,s) = ForwardDiff.derivative(s -> derivative2(l,s),s)