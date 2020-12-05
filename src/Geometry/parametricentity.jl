"""
    AbstractEntity{D,T}

A geometric entity mapping `D` into points of type `T`. 
"""
abstract type AbstractEntity{D,T} end

"""
    ParametricEntity{D,T}

A geometric entity given by an explicit `parametrization::Function` mapping `D`
into points of type `T`. 
"""
struct ParametricEntity{D,T} <: AbstractEntity{D,T}
    parametrization::Function
end

domain(p::ParametricEntity{D}) where {D<:AbstractReferenceShape} = D()

Base.eltype(p::ParametricEntity{D,T}) where {D,T} = T

ambient_dimension(p::ParametricEntity) = length(eltype(p))

geometric_dimension(p::ParametricEntity) = geometric_dimension(domain(p))

function ParametricEntity(f,d::AbstractReferenceShape)
    F = typeof(f)
    D = typeof(d)
    T = f(center(d)) |> typeof
    return ParametricEntity{D,T}(f)
end  
ParametricEntity{D}(f) where {D} = ParametricEntity(f,D())

(par::ParametricEntity)(x) = par.parametrization(x)

jacobian(psurf::ParametricEntity,s::SVector) = ForwardDiff.jacobian(psurf.parametrization,s::SVector)
jacobian(psurf::ParametricEntity,s) = jacobian(psurf,SVector(s))

function normal(ent::AbstractEntity,u)
    N = ambient_dimension(ent)
    M = geometric_dimension(ent)
    msg = "computing the normal vector requires the entity to be of co-dimension one."
    @assert N-M == 1 msg
    if M == 1 # a line in 2d
        t = jacobian(ent,u)
        ν = SVector(t[2],-t[1])
        return ν/norm(ν)
    elseif M == 2 # a surface in 3d
        j  = jacobian(ent,u)    
        t₁ = j[:,1]
        t₂ = j[:,2]
        ν  = cross(t₁,t₂)
        return ν/norm(ν)
    else
        notimplemented()    
    end
end

# TODO: move this to gmshIO.jl?
struct GmshParametricEntity{M} 
    # dim=M
    tag::Int
    model::String
    domain::HyperRectangle{M,Float64}
end

# FIXME: this constructor (or even better, an inner constructor) shoudl check
# that the gmsh surface is not trimmed.
function GmshParametricEntity(dim::Int,tag::Int,model=gmsh.model.getCurrent())
    (umin,vmin),(umax,vmax) = gmsh.model.getParametrizationBounds(dim,tag)
    rec = HyperRectangle(umin,vmin,umax-umin,vmax-vmin)
    return GmshParametricEntity{dim}(tag,model,rec,[rec])
end
GmshParametricEntity(dim::Integer,tag::Integer,args...;kwargs...) = GmshParametricEntity(Int(dim),Int(tag),args...;kwargs...)

function (par::GmshParametricEntity{N})(x) where {N}
    if N === 1
        return gmsh.model.getValue(N,par.tag,x)
    elseif N===2
        return gmsh.model.getValue(N,par.tag,[x[1],x[2]])
    else
        error("got N=$N, values must be 1 or 2")
    end
end

function jacobian(psurf::GmshParametricEntity{N},s::Point) where {N}
    if N==1
        jac = gmsh.model.getDerivative(N,psurf.tag,s)
        return reshape(jac,3,N)
    elseif N==2
        jac = gmsh.model.getDerivative(N,psurf.tag,[s[1],s[2]])
        return reshape(jac,3,N)
    else
        error("got N=$N, values must be 1 or 2")
    end
end