"""
    abstract type AbstractReferenceShape{N}
    
A reference polygon of dimension `N`.

Reference shapes are used mostly for dispatch purposes when defining instances of an 
[`AbstractElement`](@ref). 

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for 
"""
abstract type AbstractReferenceShape{N} end

"""
    struct ReferenceLine
    
Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: AbstractReferenceShape{1} 
end    
Base.in(x,::ReferenceLine) = 0 ≤ x ≤ 1

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

# TODO: generalize structs above to `RefereceSimplex{N}` and `RefereceCuboid{N}`

"""
    abstract type AbstractElement{R,N}

Abstract shape given by the image of a parametrization with domain `R::AbstractReferenceShape`.

The type parameter `N` reprenset the dimension of the embedding space. 

The dimension of the reference space can be obtained from `R`. 
"""
abstract type AbstractElement{R,N} end

"""
    abstract type PolynomialElement{R,M} <: AbstractElement{R,M}
    
Abstract elements whose parametrization is a polynomial.
"""
abstract type PolynomialElement{R,N} <: AbstractElement{R,N} end

"""
    struct LagrangeElement{R,Np,N,T}        

A lagrange elemement is reprensented as a polynomial mapping the `Np` reference lagrangian nodes of the reference element `R` into `vtx`.

The element's parametrization is fully determined by the image of the `Np` reference points, which are stored in `vtx`. 
"""
struct LagrangeElement{R,Np,N,T} <: PolynomialElement{R,N}
    vtx::SVector{Np,Point{N,T}}
end

get_order(el::LagrangeElement{ReferenceLine,Np}) where {Np} = Np + 1

# For a lagrangian element of order P, there are Np = (P+1)(P+2)/2 elements
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
        # TODO: make a macro for this in Utils?    
        error("method not (yet) implemented")
    end
end

"""
    (el::LagrangeElement)(x)

Evaluate the underlying parametrization of the element `el` at point `x`. This is the push-forward map for the element. 
"""
function (el::LagrangeElement) end

function (el::LagrangeElement{ReferenceLine,2})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine    
    vtx[1] + (vtx[2] - vtx[1]) * u    
end    

function (el::LagrangeElement{ReferenceTriangle,3})(u)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle
    vtx[1] + (vtx[2] - vtx[1]) * u[1] + (vtx[3] - vtx[1]) * u[2]
end 

function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert length(u) == 2    
    @assert u ∈ ReferenceTriangle
    vtx[1] + (vtx[2] - vtx[1]) * u[1] + (vtx[3] - vtx[1]) * u[2]
end 

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
