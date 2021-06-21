const INTERFACE_METHODS = Vector{Symbol}()

# Maybe just write the text into a file instead of using a macro here?
macro import_interface()
    ex = Expr(:block)
    ex.args = [:(import WaveProp: $f) for f in INTERFACE_METHODS]
    return ex
end

macro export_interface()
    ex = Expr(:block)
    ex.args = [:(export $f) for f in INTERFACE_METHODS]
    return ex
end

"""
    ambient_dimension(x)

Dimension of the ambient space where `x` lives. For geometrical objects this can
differ from its [`geometric_dimension`](@ref); for example a triangle in `ℝ³` has
ambient dimension `3` but geometric dimension `2`, while a curve in `ℝ³` has
ambient dimension 3 but geometric dimension 1.
"""
function ambient_dimension end
push!(INTERFACE_METHODS,:ambient_dimension)

"""
    geometric_dimension(x)
    geometric_dimension(Ω::Domain)

Number of degrees of freedom necessary to locally represent the geometrical
object. For example, lines have geometric dimension of 1 (whether in `ℝ²` or in
`ℝ³`), while surfaces have geometric dimension of 2.

When the argument is a `Domain`, return the largest geometric dimension
encoutered.
"""
function geometric_dimension end
push!(INTERFACE_METHODS,:geometric_dimension)

"""
    dimension(space)

The length of a basis for `space`; i.e. the number of linearly independent elements
required to span `space`.
"""
function dimension end
push!(INTERFACE_METHODS,:dimension)

"""
    tag(::AbstractEntity)

Integer tag commonly used to idetify geometrical entities.
"""
function tag end
push!(INTERFACE_METHODS,:tag)

"""
    boundary(ω)

Return the boundary of `ω`. For a mesh element gives the `d-1` dimensional
elements composing its boundary, while for an entity gives the corresponding
`d-1` dimensional entities.
"""
function boundary end
push!(INTERFACE_METHODS,:boundary)

"""
    diameter(Ω)

Largest distance between `x` and `y` for `x,y ∈ Ω`.
"""
function diameter end
push!(INTERFACE_METHODS,:diameter)

"""
    radius(Ω)

Half the [`diameter`](@ref).
"""
function radius end
push!(INTERFACE_METHODS,:radius)

"""
    center(Ω)
"""
function center end
push!(INTERFACE_METHODS,:center)

"""
    return_type(f)

The type returned by the function-like object `f`.
"""
function return_type end
push!(INTERFACE_METHODS,:return_type)

"""
    jacobian(f,x)

Given a (possibly vector-valued) function `f : 𝐑ᵐ → 𝐅ᵐ`, return the `m × n`
matrix `Aᵢⱼ = ∂fᵢ/∂x̂ⱼ`.
"""
function jacobian end
push!(INTERFACE_METHODS,:jacobian)

"""
    normal(el,x̂)
    normal(jac::SMatrix)

The unit normal vector at coordinate `x̂`, guaranteed to be orthogonal to all
columns of `jacobian(el,x)`.
"""
function normal end
push!(INTERFACE_METHODS,:normal)

"""
    domain(f)

The domain of the function `f`, meaning
"""
function domain end
push!(INTERFACE_METHODS,:domain)

"""
    parametrization(el)

Return the underlying parametrization of `el`.
"""
function parametrization end
push!(INTERFACE_METHODS,:parametrization)

"""
    entities(Ω::Domain)
    entities(M::AbstractMesh)

Return the geometrical entities composing `Ω`.
"""
function entities end
push!(INTERFACE_METHODS,:entities)

"""
    bounding_box(data)

Create an axis-aligned bounding box containing all of `data`.
"""
function bounding_box end
push!(INTERFACE_METHODS,:bounding_box)


"""
    reference_nodes(::LagrangeElement{D,Np,T})

Return the reference nodes on `D` used for the polynomial interpolation. The
function values on these nodes completely determines the interpolating
polynomial.

We use the same convention as `gmsh` for defining the reference nodes and their
order (see [node
ordering](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) on `gmsh`
documentation).
"""
function reference_nodes end
push!(INTERFACE_METHODS,:reference_nodes)
