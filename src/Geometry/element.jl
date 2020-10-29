"""
    abstract type AbstractReferenceShape{N}
    
A reference polygon of dimension `N`.

Reference shapes are used mostly for dispatch purposes when defining
instances of an [`AbstractElement`](@ref).

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for
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

"""
    struct ReferenceTetrahedron
    
Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: AbstractReferenceShape{3} 
end    
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1] && 0 ≤ x[3] ≤ 1 - x[1] - x[2]

# TODO: generalize structs above to `ReferenceSimplex{N}` and `ReferenceCuboid{N}`

"""
    abstract type AbstractElement{R,N}

Abstract shape given by the image of a parametrization with domain
`R<:AbstractReferenceShape`.

The type parameter `N` represents the dimension of the embedding space.

The dimension of the reference space can be obtained from `R`. 
"""
abstract type AbstractElement{R,N} end

"""
    reference_element(el::AbstractElement)

Return an instance of the singleton type `R`; i.e. the reference element.
"""
domain(::Type{<:AbstractElement{R}}) where {R} = R()
domain(el::AbstractElement) = domain(typeof(el))

function normal(el::AbstractElement,u)
    N = ambient_dimension(el)
    M = geometric_dimension(el)
    msg = "computing the normal vector requires the element to be of co-dimension one."
    @assert N-M == 1 msg
    if M == 1 # a line in 2d
        t⃗ = jacobian(el,u)
        n⃗ = SVector(-t⃗[2],t⃗[1])
        return n⃗/norm(n⃗)
    elseif M == 2 # a surface in 3d
        j  = jacobian(el,u)    
        t⃗₁ = j[:,1]
        t⃗₂ = j[:,2]
        n⃗  = cross(t⃗₁,t⃗₂)
        return n⃗/norm(n⃗)
    else
        @notimplemented    
    end            
end    

push_forward_map(el,x̂) = map(el,x̂)

function push_forward_map(el,x̂,ŵ)
    x   = push_forward_map(el,x̂)
    w   = map(zip(x̂,ŵ)) do (x̂,ŵ)
        jac = jacobian(el,x̂)
        g   = transpose(jac)*jac |> det
        sqrt(g)*prod(ŵ)
    end 
    return x,w
end   
# FIXME: the function above is somewhat inefficient when the ambient and
# geometric dimensions of the element are the same. In that case `μ` simplifies
# to the usual `|det(jac)|`. This should be easy to fix by checking e.g. whether
# `jac` is a square matrix. Since these are static arrays there should be no
# runtime overhead compared to the hand-written version

"""
    geometric_dimension(el::AbstractElement{R})

Return the geometric dimension of `el`, i.e. the number of variables needed to locally
parametrize the element.
"""
geometric_dimension(t::Type{<:AbstractElement}) = geometric_dimension(domain(t))
geometric_dimension(el)                         = geometric_dimension(typeof(el))

"""
    ambient_dimension(el::AbstractElement)

Return the dimension of the ambient space where `el` lives.
"""
ambient_dimension(el::AbstractElement{R,N}) where {R,N} = N
ambient_dimension(::Type{<:AbstractElement{R,N}}) where {R,N} = N

"""
    abstract type PolynomialElement{R,M} <: AbstractElement{R,M}
    
Abstract elements whose parametrization is a polynomial.
"""
abstract type PolynomialElement{R,N} <: AbstractElement{R,N} end

"""
    struct LagrangeElement{R,Np,N,T}        
    
Fields:
    - `nodes`

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

get_nodes(el::LagrangeElement) = el.nodes

get_order(el::LagrangeElement{ReferenceLine,Np}) where {Np} = Np + 1

dim(el::LagrangeElement{R,Np,N}) where {R,Np,N} = N

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
        @notimplemented
    end
end

"""
    (el::LagrangeElement)(x)

Evaluate the underlying parametrization of the element `el` at point `x`. This is the push-forward map for the element. 
"""
function (el::LagrangeElement) end

function (el::LagrangeElement{ReferenceLine,2})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    nodes = get_nodes(el)
    nodes[1] + (nodes[2] - nodes[1]) * u[1]    
end    

function (el::LagrangeElement{ReferenceTriangle,3})(u)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle()
    nodes = get_nodes(el)
    nodes[1] + (nodes[2] - nodes[1]) * u[1] + (nodes[3] - nodes[1]) * u[2]
end 

function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert length(u) == 3    
    @assert u ∈ ReferenceTetrahedron()
    nodes = get_nodes(el)
    nodes[1] + (nodes[2] - nodes[1]) * u[1] + (nodes[3] - nodes[1]) * u[2] + (nodes[4] - nodes[1]) * u[3]
end 

"""
    gradient(el::LagrangeElement,x)

Evaluate the gradient of the underlying parametrization of the element `el` at point `x`. 
"""
function jacobian(el::LagrangeElement{ReferenceLine,2}, u)
    N = dim(el)    
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()    
    nodes = get_nodes(el)
    return SMatrix{N,1}((nodes[2] - nodes[1])...)
end    

function jacobian(el::LagrangeElement{ReferenceTriangle,3}, u) 
    N = dim(el)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle()
    nodes = get_nodes(el)
    SMatrix{N,2}(
        (nodes[2] - nodes[1])..., 
        (nodes[3] - nodes[1])...
    )
end 

function jacobian(el::LagrangeElement{ReferenceTetrahedron,4}, u)
    N = dim(el)    
    @assert length(u) == 3   
    @assert u ∈ ReferenceTriangle()
    nodes = get_nodes(el)
    SMatrix{N,3}( 
        (nodes[2] - nodes[1])...,
        (nodes[3] - nodes[1])...,
        (nodes[4] - nodes[1])...
    )
end 

# define some aliases for convenience
const LagrangeLine        = LagrangeElement{ReferenceLine}
const LagrangeTriangle    = LagrangeElement{ReferenceTriangle}
const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}

"""
    abstract type ParametricElement{R,N} <: AbstractElement{R,N}
    
An element given by some explicit parametric function.        
"""
abstract type AbstractParametricElement{R,N} <: AbstractElement{R,N} end

"""
    ParametricElement{R,N,F}

An element represented as the explicit mapping `f::F` with domain `R`.  
"""
struct ParametricElement{R,N,F} <: AbstractParametricElement{R,N}
    parametrization::F    
end

(el::ParametricElement)(u) = el.parametrization(u)

