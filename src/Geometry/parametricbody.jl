"""
    abstract type AbstractParametricBody <: AbstractEntity
    
An `AbstractEntity` for which a parametric representation of the boundary is
available. 
"""
abstract type AbstractParametricBody <: AbstractEntity end

boundary(ent::AbstractParametricBody) = ent.boundary

function ambient_dimension(ent::AbstractParametricBody)
    bnd = boundary(ent)    
    d   = ambient_dimension(first(bnd))
    msg = "all entities on the boundary of a body must 
            have the same ambient dimension"
    @assert all(b->ambient_dimension(b)==d,bnd) msg
    return d
end    

"""
    ClosedEntity <: AbstractParametricBody

A geometric entity whose boundary is defined parametrically.
"""
struct ClosedEntity <: AbstractParametricBody
    dim::UInt8
    tag::Int
    boundary::Vector{ParametricEntity}
    function ClosedEntity(d, t, bnd)
        ent = new(d, t, bnd)
        _global_add_entity!(ent)
        return ent
    end        
end

function ClosedEntity(boundary::Vector{<:ParametricEntity})
    dims = geometric_dimension.(boundary) 
    d    = dims[1]
    @assert all(i->i==d,dims)
    d = d+1 # dimension of body is one greater than dimension of its boundary
    t = _new_tag(d)
    ClosedEntity(d, t, boundary)
end    
ClosedEntity(boundary::ParametricEntity) = ClosedEntity([boundary,])
ClosedEntity(;boundary::Vector{<:ParametricEntity}) = ClosedEntity(boundary)


key(ent::ClosedEntity) = (ent.dim,ent.tag)
geometric_dimension(ent::ClosedEntity) = ent.dim

struct Circle <: AbstractParametricBody
    # dim = 2
    tag::Int
    center::SVector{2,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Circle(t,c,r,bnd)
        d = 2
        ent = new(t,c,r,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Circle(;center=(0, 0),radius=1)
    f          = (s) -> center .+ radius .* SVector(cospi(2 * s[1]), sinpi(2 * s[1]))
    domain     = ReferenceLine()
    ent        = ParametricEntity(f, domain)
    tag        = _new_tag(2) # generate a unique tag for entities of dimension 2
    return Circle(tag,center, radius, [ent])
end
Base.in(pt,circ::Circle) = norm(pt .- circ.center) < circ.radius
key(ent::Circle) = (2,ent.tag)
geometric_dimension(ent::Circle) = 2

struct Kite <: AbstractParametricBody
    # dim = 2
    tag::Int    
    center::SVector{2,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Kite(t,c,r,bnd)
        d   = 2
        ent = new(t,c,r,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Kite(;radius=1,center=(0, 0))
    f = (s) -> center .+ radius .* SVector(cospi(2 * s[1]) + 0.65 * cospi(4 * s[1]) - 0.65,
                              1.5 * sinpi(2 * s[1]))
    domain = ReferenceLine()
    surf   = ParametricEntity(f, domain)
    t      = _new_tag(2)
    return Kite(t,center, radius, [surf])
end
key(ent::Kite) = (2,ent.tag)
geometric_dimension(ent::Kite) = 2

struct Droplet <: AbstractParametricBody
    # dim = 2
    tag::Int    
    center::SVector{2,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Droplet(t,c,r,bnd)
        d   = 2
        ent = new(t,c,r,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Droplet(;radius=1,center=(0,0))
    f = (s) -> center .+ radius .* SVector(2 * sinpi(s[1]),-sinpi(2*s[1]))
    domain = ReferenceLine()
    surf   = ParametricEntity(f, domain)
    t      = _new_tag(2)
    return Droplet(t,center, radius, [surf])
end    
key(ent::Droplet) = (2,ent.tag)
geometric_dimension(ent::Droplet) = 2

struct Boomerang <: AbstractParametricBody
    # dim = 2
    tag::Int    
    center::SVector{2,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Boomerang(t,c,r,bnd)
        d   = 2
        ent = new(t,c,r,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Boomerang(;radius=1,center=(0,0))
    f = (s) -> center .+ radius .* SVector(-2/3 * sinpi(3*s[1]),-sinpi(2*s[1]))
    domain = ReferenceLine()
    surf   = ParametricEntity(f, domain)
    t      = _new_tag(2)
    return Droplet(t,center, radius, [surf])
end    
key(ent::Boomerang) = (2,ent.tag)
geometric_dimension(ent::Boomerang) = 2

# struct Ellipsoid{T} <: AbstractParametricBody{3,2,T}
#     center::SVector{3,T}
#     paxis::SVector{3,T}
#     parts::Vector{ParametricEntity{3,2,T}}
# end

# function Ellipsoid{T}(;center=zeros(3),paxis=ones(3)) where {T}
#     nparts = 6
#     domain = HyperRectangle((-1.,-1.),(1.,1.))
#     parts  = Vector{ParametricEntity}(undef,nparts)
#     for id=1:nparts
#         param(x)     = _ellipsoid_parametrization(x[1],x[2],id,paxis,center)
#         parts[id]    = ParametricEntity(param,domain,[domain])
#     end
#     return Ellipsoid{T}(center,paxis,parts)
# end
# Ellipsoid(args...;kwargs...) = Ellipsoid{Float64}(args...;kwargs...)

struct Sphere <: AbstractParametricBody
    # dim = 3
    tag::Int
    center::SVector{3,Float64}
    radius::Float64
    boundary::Vector{ParametricEntity}
    function Sphere(t,c,r,bnd)
        d = 3
        ent = new(t,c,r,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Sphere(;center=(0, 0,0),radius=1)
    nparts = 6
    domain = ReferenceSquare()
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param     = (x) -> _sphere_parametrization(x[1],x[2],id,radius,center)
        parts[id] = ParametricEntity(param,domain)
    end
    tag = _new_tag(3)
    return Sphere(tag,center,radius,parts)
end
Base.in(pt,sph::Sphere) = norm(pt .- sph.center) < sph.radius
key(ent::Sphere) = (3,ent.tag)
geometric_dimension(ent::Sphere) = 3

struct Cube <: AbstractParametricBody
    # dim = 3
    tag::Int
    origin::SVector{3,Float64}
    paxis::SVector{3,Float64}
    boundary::Vector{ParametricEntity}
    function Cube(t,c,p,bnd)
        d = 3
        ent = new(t,c,p,bnd)    
        _global_add_entity!(ent)
        return ent
    end    
end

function Cube(;origin=SVector(0,0,0),paxis=(1,1,1))
    nparts = 6
    domain = ReferenceSquare()
    parts  = Vector{ParametricEntity}(undef,nparts)
    for id=1:nparts
        param     = (x) -> _cube_parametrization(x[1],x[2],id,paxis,origin)
        parts[id] = ParametricEntity(param,domain)
    end
    tag = _new_tag(3)
    return Cube(tag,origin,paxis,parts)
end
# Base.in(pt,cub::Cube) = norm(pt .- sph.center) < sph.radius
key(ent::Cube) = (3,ent.tag)
geometric_dimension(ent::Cube) = 3

# struct Bean{T} <: AbstractParametricBody{3,2,T}
#     center::SVector{3,T}
#     paxis::SVector{3,T}
#     parts::Vector{ParametricEntity{3,2,T}}
# end
# function Bean{T}(;center=zeros(3),paxis=ones(3)) where {T}
#     nparts = 6
#     domain = HyperRectangle(-1.,-1.,2.,2.)
#     parts  = Vector{ParametricEntity}(undef,nparts)
#     for id=1:nparts
#         param     = (x) -> _bean_parametrization(x[1],x[2],id,paxis,center)
#         parts[id] = ParametricEntity(param,domain,[domain])
#     end
#     return Bean{T}(center,paxis,parts)
# end
# Bean(args...;kwargs...) = Bean{Float64}(args...;kwargs...)

# struct Acorn{T} <: AbstractParametricBody{3,2,T}
#     center::SVector{3,T}
#     radius::T
#     rotation::SVector{3,T}
#     parts::Vector{ParametricEntity{3,2,T}}
# end
# function Acorn{T}(;center=zeros(3),radius=1.0,rotation=(0,0,0)) where {T}
#     nparts = 6
#     domain = HyperRectangle(-1.,-1.,2.,2.)
#     parts  = Vector{ParametricEntity}(undef,nparts)
#     for id=1:nparts
#         param        = (x) -> _acorn_parametrization(x[1],x[2],id,radius,center,rotation)
#         parts[id]    = ParametricEntity(param,domain,[domain])
#     end
#     return Acorn{T}(center,radius,rotation,parts)
# end
# Acorn(args...;kwargs...) = Acorn{Float64}(args...;kwargs...)

# struct Cushion{T} <: AbstractParametricBody{3,2,T}
#     center::SVector{3,T}
#     radius::T
#     rotation::SVector{3,T}
#     parts::Vector{ParametricEntity{3,2,T}}
# end
# function Cushion{T}(;center=zeros(3),radius=1.0,rotation=(0,0,0)) where {T}
#     nparts = 6
#     domain = HyperRectangle(-1.,-1.,2.,2.)
#     parts  = Vector{ParametricEntity}(undef,nparts)
#     for id=1:nparts
#         param     = (x) -> _cushion_parametrization(x[1],x[2],id,radius,center,rotation)
#         parts[id] = ParametricEntity(param,domain,[domain])
#     end
#     return Cushion{T}(center,radius,rotation,parts)
# end
# Cushion(args...;kwargs...) = Cushion{Float64}(args...;kwargs...)

function _cube_parametrization(u,v,id,paxis=SVector(1,1,1),origin=SVector(0.0,0.0,0.0))
    if id==1
        x = SVector(1.,u,v)
    elseif id==2
        x = SVector(u,1.,v)
    elseif id==3
        x = SVector(u,v,1.)
    elseif id==4
        x = SVector(0.,u,v)
    elseif id==5
        x = SVector(u,0.,v)
    elseif id==6
        x = SVector(u,v,0.)
    end
    return origin + paxis .* x
end

function _sphere_parametrization(u,v,id,rad=1,center=zeros(3))
    # map [0,1]×[0,1] to [-1,1] × [-1,1]
    u = 2u-1
    v = 2v-1
    # parametrization of 6 patches
    if id==1
        x = SVector(1.,u,v)
    elseif id==2
        x = SVector(-u,1.,v)
    elseif id==3
        x = SVector(u,v,1.)
    elseif id==4
        x = SVector(-1.,-u,v)
    elseif id==5
        x = SVector(u,-1.,v)
    elseif id==6
        x = SVector(-u,v,-1.)
    end
    return center .+ rad.*x./sqrt(u^2+v^2+1)
end

# function _ellipsoid_parametrization(u,v,id,paxis=ones(3),center=zeros(3))
#     x = _sphere_parametrization(u,v,id)
#     return x .* paxis .+ center
# end

# function _bean_parametrization(u,v,id,paxis=one(3),center=zeros(3))
#     x = _sphere_parametrization(u,v,id)
#     a = 0.8; b = 0.8; alpha1 = 0.3; alpha2 = 0.4; alpha3=0.1
#     x[1] = a*sqrt(1.0-alpha3*cospi(x[3])).*x[1]
#     x[2] =-alpha1*cospi(x[3])+b*sqrt(1.0-alpha2*cospi(x[3])).*x[2]
#     x[3] = x[3];
#     return x .* paxis .+ center
# end

# function _acorn_parametrization(u,v,id,radius,center,rot)
#     Rx = [1 0 0;0 cos(rot[1]) sin(rot[1]);0 -sin(rot[1]) cos(rot[1])]
#     Ry = [cos(rot[2]) 0 -sin(rot[2]);0 1 0;sin(rot[2]) 0 cos(rot[2])]
#     Rz = [cos(rot[3]) sin(rot[3]) 0;-sin(rot[3]) cos(rot[3]) 0;0 0 1]
#     R  = Rz*Ry*Rx;
#     x = _sphere_parametrization(u,v,id)
#     th,phi,_ = cart2sph(x...)
#     r=0.6+sqrt(4.25+2*cos(3*(phi+pi/2)))
#     x[1] = r.*cos(th).*cos(phi)
#     x[2] = r.*sin(th).*cos(phi)
#     x[3] = r.*sin(phi)
#     x    = R*x
#     return radius.*x .+ center
# end

# function _cushion_parametrization(u,v,id,radius,center,rot)
#     Rx = [1 0 0;0 cos(rot[1]) sin(rot[1]);0 -sin(rot[1]) cos(rot[1])]
#     Ry = [cos(rot[2]) 0 -sin(rot[2]);0 1 0;sin(rot[2]) 0 cos(rot[2])]
#     Rz = [cos(rot[3]) sin(rot[3]) 0;-sin(rot[3]) cos(rot[3]) 0;0 0 1]
#     R  = Rz*Ry*Rx;
#     x = _sphere_parametrization(u,v,id)
#     th,phi,_ = cart2sph(x...)
#     r = sqrt(0.8+0.5*(cos(2*th)-1).*(cos(4*phi)-1));
#     x[1] = r.*cos(th).*cos(phi)
#     x[2] = r.*sin(th).*cos(phi)
#     x[3] = r.*sin(phi)
#     x    = R*x
#     return radius.*x .+ center
# end

# function cart2sph(x,y,z)
#     azimuth = atan(y,x)
#     elevation = atan(z,sqrt(x^2 + y^2))
#     r = sqrt(x^2 + y^2 + z^2)
#     return azimuth, elevation, r
# end


# struct GmshParametricBody{M} <: AbstractParametricBody{3,M,Float64}
#     parts::Vector{GmshParametricEntity{M}}
# end

# function GmshParametricBody(dim,tag,model=gmsh.model.getCurrent())
#     dimtags = gmsh.model.getBoundary((dim,tag))
#     body    = GmshParametricBody{2}([])
#     for dimtag in dimtags
#         patch = GmshParametricEntity(Int(dimtag[1]),Int(dimtag[2]),model)
#         push!(body.parts,patch)
#     end
#     return body
# end
