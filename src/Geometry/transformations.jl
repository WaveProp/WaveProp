"""
    abstract type GeometricTransformation
    
Abstract type representing a change of variable mapping `domain` to `range`

In order to aid with various methods, instances `transf::GeometricTransformation` are expected to implement
- `domain`
- `range`
- `(::GeometricTransformation)(u)`
- `jacobian(::GeometricTransformation,u)`
"""
abstract type GeometricTransformation end

"""
    struct IMT{A,P} <: GeometricTransformation
    
One-dimensional change of variables mapping `[0,1] -> [0,1]` with the property that 
all derivatives of the jacobian vanish at `0`.

# FIXME: how do you add a reference to a docstring?
See [Davis and Rabinowitz](https://www.elsevier.com/books/methods-of-numerical-integration/davis/978-0-12-206360-2)
"""
struct IMT{A,P} <: GeometricTransformation
end

domain(::IMT) = ReferenceLine()
range(::IMT)  = ReferenceLine()

IMT(;a=1,p=1) = IMT{a,p}()

@fastmath function (f::IMT{A,P})(x) where {A,P}
    @. exp(A*(1-1/x^P))
end

@fastmath derivative(f::IMT{A,P},x) where {A,P} = @. f(x) * A*P*1 / x^(P+1)
jacobian(f::IMT,x) = derivative(f,x) |> SMatrix{1,1}

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
struct Duffy{N} <: GeometricTransformation end

domain(::Duffy{2}) = ReferenceSquare()
range(::Duffy{2})  = ReferenceTriangle()

function (::Duffy{2})(u)
    SVector(u[1],(1-u[1])*u[2])
end    

function jacobian(::Duffy{2},u)
    SMatrix{2,2,Float64}(1,0,-u[2],(1-u[1]))
end    