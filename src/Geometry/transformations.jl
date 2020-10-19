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

"""
    struct Duffy{N} <: ChangeOfVariables
    
Change of variables mapping the reference `HyperCube{N}` to the `Simplex{N}`
with the property that the jacobian vanishes at the `𝐞₁` vertex of the simplex.
"""
struct Duffy{N} <: GeometricTransformation end

domain(::Duffy{2}) = ReferenceSquare()
range(::Duffy{2})  = ReferenceTriangle()

function (::Duffy{2})(u::Point{2})
    SVector(u[1],(1-u[1])*u[2])
end    

function jacobian(::Duffy{2},u::Point{2})
    SMatrix{2,2,Float64}(1,0,-u[2],(1-u[1]))
end    