"""
    abstract type AbstractEntity
    
Entity of geometrical nature. Identifiable throught its `(dim,tag)` key.
"""
abstract type AbstractEntity end

"""
    struct ElementaryEntity <: AbstractEntity
    
The most basic representation of an `AbstractEntity`. 

# Fields:
- `dim::UInt8`: the geometrical dimension of the entity (e.g. line has `dim=1`, surface has `dim=2`, etc)
- `tag::Int64`: an integer tag associated to the entity
- `boundary::Vector{ElementaryEntity}`: the entities of dimension `dim-1` forming the entity's boundary  
"""
struct ElementaryEntity <: AbstractEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{ElementaryEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{ElementaryEntity})
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d-1 msg
        end
        ent = new(d, t, boundary)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities. 
        _global_add_entity!(ent)
        return ent
    end
end

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim,tag)
    ElementaryEntity(dim,tag,ElementaryEntity[])
end
ElementaryEntity(dim) = ElementaryEntity(dim,_new_tag(dim))

geometric_dimension(ω::ElementaryEntity) = ω.dim

"""
    tag(ω::ElementaryEntity)

Return the unique tag (for a given dimension) of the elementary entity.
"""
key(ω::ElementaryEntity) = (geometric_dimension(ω), ω.tag)

"""
    boundary(ω::ElementaryEntity)

Return the vector of elementary entities making the boundary.
"""
boundary(ω::ElementaryEntity) = ω.boundary

"""
    ==(Ω1::ElementaryEntity,Ω2::ElementaryEntity) 

Two elementary entities are considered equal if their `dim` and `tag` fields
match.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
define it, and therefore global variables like [`TAGS`](@ref) are needed to make
sure newly created `AbstractEntities` have a new `(dim,tag)` identifier. 
"""
function Base.:(==)(Ω1::AbstractEntity, Ω2::AbstractEntity)
    d1,t1 = key(Ω1)
    d2,t2 = key(Ω2)
    d1 == d2  || return false
    abs(t1) == abs(t2) || return false
    # boundary(Ω1) == boundary(Ω2) || return false # this should not be needed
    return true
end


"""
    ParametricEntity{D,F}

A geometric entity given by an explicit `parametrization::F` with domain given
by a singleton type `D`.

This differs from an `ElementaryEntity` in that the underlying representation of
the geometric entity is available. 
"""
struct ParametricEntity{D,F} <: AbstractEntity
    dim::UInt8
    tag::Int 
    parametrization::F
    function ParametricEntity{D,F}(dim,tag,f::F) where {D,F}
        ent = new{D,F}(dim,tag,f)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities. 
        _global_add_entity!(ent)        
        return ent
    end    
end

"""
    const ParametricCurve{F} = ParametricEntity{ReferenceLine,F}

A curve defined by the function `f::F`.
"""
const ParametricCurve{F} = ParametricEntity{ReferenceLine,F} 

function Base.reverse(ent::ParametricEntity{D,F}) where {D,F}
    ParametricEntity{D,F}(ent.dim,-ent.tag,ent.parametrization)
end

function ParametricEntity{D,F}(f::F) where {D,F}
    d = geometric_dimension(D)
    t = _new_tag(d) # automatically generate a new (valid) tag
    ParametricEntity{D,F}(d,t,f)
end       

function ParametricEntity(f,d::AbstractReferenceShape)
    F = typeof(f)
    D = typeof(d)
    return ParametricEntity{D,F}(f)
end  
ParametricEntity{D}(f) where {D} = ParametricEntity(f,D())

domain(p::ParametricEntity{D}) where {D<:AbstractReferenceShape} = D()
parametrization(p::ParametricEntity) = p.parametrization 

key(ent::ParametricEntity) = (ent.dim,ent.tag)

# FIXME: what should this be?
boundary(ent::ParametricEntity) = ElementaryEntity[]

function Base.eltype(p::ParametricEntity)
    # NOTE: this relies on the brittle promote_op
    d = domain(p)    
    x = center(d)
    f = parametrization(p)
    T = Base.promote_op(f,typeof(x))
end    

ambient_dimension(p::ParametricEntity)   = length(eltype(p))
geometric_dimension(p::ParametricEntity) = geometric_dimension(domain(p))

(par::ParametricEntity)(x) = par.parametrization(x)

jacobian(psurf::ParametricEntity,s::SVector)  = ForwardDiff.jacobian(psurf.parametrization,s::SVector)
jacobian(psurf::ParametricEntity,s)           = jacobian(psurf,SVector(s))

function normal(ent::AbstractEntity,u)  
    dim,tag = key(ent)
    s = sign(tag)
    N    = ambient_dimension(ent)
    M    = geometric_dimension(ent)
    msg  = "computing the normal vector requires the entity to be of co-dimension one."
    @assert N-M == 1 msg
    if M == 1 # a line in 2d
        t = jacobian(ent,u)
        ν = SVector(t[2],-t[1])
        return s*ν/norm(ν)
    elseif M == 2 # a surface in 3d
        j  = jacobian(ent,u)    
        t₁ = j[:,1]
        t₂ = j[:,2]
        ν  = s*cross(t₁,t₂)
        return ν/norm(ν)
    else
        notimplemented()    
    end
end

"""
    line(a,b)

Create a straight line connecting points `a` and `b`. This returns an instance
of [`ParametricCurve`](@ref).
"""
function line(a::SVector,b::SVector) 
    f = (u) -> a + u[1]*(b-a)
    ParametricCurve(f)
end
line(a,b) = line(SVector(a),SVector(b))

"""
    Point{N,T} <: AbstractEntity
"""
struct Point{N,T} <: AbstractEntity
    coords::SVector{N,T}
end    

geometric_dimension(::Type{<:Point}) = 0
geometric_dimension(pt::Point) = geometric_dimension(typeof(pt))
ambient_dimension(::Point{N}) where N = N
ambient_dimension(::Type{<:Point{N}}) where N = N

geometric_dimension(::Type{<:SVector}) = 0
geometric_dimension(pt::SVector) = geometric_dimension(typeof(pt))
ambient_dimension(::SVector{N}) where N = N
ambient_dimension(::Type{<:SVector{N}}) where N = N

#####################################################################

# Variables and functions to globally keep track of entities

#####################################################################

"""
    const TAGS::OrderedDict{Int,Vector{Int}}

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""    
const TAGS     = OrderedDict{Int,Vector{Int}}()

"""
    const ENTITIES

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""    
const ENTITIES = OrderedDict{Tuple{Int,Int},AbstractEntity}()

function _global_add_entity!(ent::AbstractEntity)
    d,t = key(ent)
    _add_tag!(d,t) # add this tag to global list to make sure it is not used again
    msg = "overwriting ENTITIES: value in key ($d,$t) will be replaced"
    haskey(ENTITIES,(d,t)) && (@warn msg)
    ENTITIES[(d,t)] = ent
    return d,t
end    

"""
    _new_tag(dim)

Generate a unique tag for an `AbstractEntity` of dimension `dim`.

The implementation consists of adding one to the maximum value of `TAGS[dim]`

# See also: [`TAGS`](@ref).

"""
function _new_tag(dim)
    if !haskey(TAGS,dim) 
        return 1
    else
        tnew = maximum(TAGS[dim]) + 1
        return tnew    
    end
end

function _add_tag!(dim,tag)
    if is_new_tag(dim,tag)
        # now add key
        if haskey(TAGS,dim)   
            push!(TAGS[dim],tag)
        else
            TAGS[dim] = [tag,]
        end
    else
        # print warning but don't add duplicate tag
        msg  = "entity of dimension $dim and tag $tag already exists in TAGS. 
        Creating a possibly duplicate entity."        
        @warn msg
    end    
    return TAGS    
end 

function is_new_tag(dim,tag)
    if haskey(TAGS,dim)
        existing_tags = TAGS[dim]
        if in(tag,existing_tags) 
            return false
        end
    end
    return true
end  

clear_tags!() = empty!(TAGS)
clear_entities!() = empty!(ENTITIES)
function clear!()
    clear_tags!()
    clear_entities!()
    nothing
end

