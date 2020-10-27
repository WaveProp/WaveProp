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
    reverse(el::LagrangeLine)

Create a new element with the orientation reversed.

Because of the implicit ordering of the elements, where the boundary nodes go
first, followed by the internal nodes, the `reverse` function is slighly more
complicated than just reversing the order of the nodes.
"""
function Base.reverse(el::LagrangeLine)
    E = typeof(el)    
    x = nodes(el)
    N = length(x)
    perm = svector(N) do i
        i == 1  && return 2
        i == 2  && return 1
        return N + 3 - i
    end
    return E(x[perm])
end    

"""
    translate(el::LagrangeLine,d⃗)

Create a new element correspoding to a translation of `el` along `d⃗`
"""
function translate(el::LagrangeElement,d⃗)
    E = typeof(el)    
    x = map(x->x+d⃗,nodes(el))    
    return E(x)
end    
