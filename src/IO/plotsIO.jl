"""
    struct PlotPoints

Structure used for dispatching `SVector` to plot recipes without type-piracy.
"""
struct PlotPoints end

@recipe function f(::PlotPoints, pts::SVector{N1, SVector{N2, Float64}}) where {N1, N2}
    if N2 == 2
        xx = SVector{N1, Float64}(pt[1] for pt in pts)
        yy = SVector{N1, Float64}(pt[2] for pt in pts)
        return xx,yy
    elseif N2 == 3
        xx = SVector{N1, Float64}(pt[1] for pt in pts)
        yy = SVector{N1, Float64}(pt[2] for pt in pts)
        yy = SVector{N1, Float64}(pt[3] for pt in pts)
        return xx,yy,zz
    else
        notimplemented()
    end
end

@recipe function f(::PlotPoints,pts::AbstractVector{<:SVector{N}}) where {N}
    if N == 2
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        return xx,yy
    elseif N == 3
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        zz = [pt[3] for pt in pts]
        return xx,yy,zz
    else
        notimplemented()
    end
end

@recipe function f(::PlotPoints,pts::AbstractMatrix{<:SVector})
    PlotPoints(),vec(pts)
end

# plot a hyperrectangle
@recipe function f(rec::HyperRectangle{N}) where {N}
    seriestype := :path
    linecolor --> :black
    linestyle --> :solid
    label --> ""
    if N == 2
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1)
    elseif N == 3
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        z1, z2 = pt1[3], pt2[3]
        # upper and lower faces
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1), SVector(z1,z1,z1,z1,z1)
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1), SVector(z2,z2,z2,z2,z2)
        # lines connecting faces
        @series SVector(x1,x1), SVector(y1,y1), SVector(z1,z2)
        @series SVector(x2,x2), SVector(y1,y1), SVector(z1,z2)
        @series SVector(x2,x2), SVector(y2,y2), SVector(z1,z2)
        @series SVector(x1,x1), SVector(y2,y2), SVector(z1,z2)
    end
end

# plot domain --> plot all of its entities
@recipe function f(Ω::Domain)
    for ent in entities(Ω)
        @series begin
            ent
        end
    end
end

# recipe for parametric line
@recipe function f(ent::ParametricEntity;h=0.1)
    d =  domain(ent)
    lc = low_corner(d)
    hc = high_corner(d)
    if d isa HyperRectangle{1}
        par = ent.parametrization
        grid   --> false
        aspect_ratio --> :equal
        s       = lc[1]:h:hc[1]
        pts     = [par(v) for v in s]
        x       = [pt[1] for pt in pts]
        y       = [pt[2] for pt in pts]
        x,y
    elseif d isa HyperRectangle{2}
        legend --> false
        grid   --> false
        # aspect_ratio --> :equal
        seriestype := :surface
        xrange = lc[1]:h:hc[1]
        yrange = lc[2]:h:hc[2]
        pts    = [ent.parametrization((x,y)) for x in xrange, y in yrange]
        x      =  [pt[1] for pt in pts]
        y      =  [pt[2] for pt in pts]
        z      =  [pt[3] for pt in pts]
        x,y,z
    end
end

@recipe function f(ents::Vector{ParametricEntity};h=0.1)
    grid   --> false
    aspect_ratio --> :equal
    for ent in ents
        @series begin
            ent
        end
    end
end

# # recipe for paramatric surface
# @recipe function f(ent::ParametricEntity{ReferenceSquare};h=0.1)
#     legend --> false
#     grid   --> false
#     # aspect_ratio --> :equal
#     seriestype := :surface
#     xrange = 0:h:1
#     yrange = 0:h:1
#     pts    = [ent.parametrization((x,y)) for x in xrange, y in yrange]
#     x      =  [pt[1] for pt in pts]
#     y      =  [pt[2] for pt in pts]
#     z      =  [pt[3] for pt in pts]
#     x,y,z
# end

# recipe for parametric body
# @recipe function f(bdy::AbstractParametricBody)
#     label --> ""
#     # aspect_ratio --> :equal
#     for patch in boundary(bdy)
#         @series begin
#             patch
#         end
#     end
# end

# # recipe for many parametric bodies
# @recipe function f(bdies::Vector{<:AbstractParametricBody})
#     aspect_ratio --> :equal
#     for bdy in bdies
#         @series begin
#             bdy
#         end
#     end
# end

@recipe function f(mesh::NystromMesh)
    label --> ""
    grid   --> false
    aspect_ratio --> :equal
    for (E,els) in mesh.elements
        for el in els
            @series begin
                el
            end
        end
    end
end

@recipe function f(mesh::GenericMesh,Ω::Domain)
    view(mesh,Ω)
end

@recipe function f(mesh::SubMesh)
    label --> ""
    grid   --> false
    aspect_ratio --> :equal
    for E in keys(mesh)
        for el in ElementIterator(mesh,E)
            @series begin
                el
            end
        end
    end
end

# FIXME: write better recipes
@recipe function f(el::LagrangeTriangle)
    label --> ""
    vtx = nodes(el)
    for n in 1:3
        is = n
        ie = 1 + (n%3)
        @series begin
            [vtx[is],vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeSquare)
    label --> ""
    vtx = nodes(el)
    for n in 1:4
        is = n
        ie = 1 + (n%4)
        @series begin
            [vtx[is],vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeLine)
    vtx = nodes(el)
    [vtx[1],vtx[2]]
end

@recipe function f(el::AbstractElement{ReferenceLine};npts=10)
    label --> ""
    grid   --> false
    s       =  LinRange(0,1,npts)
    pts     = [el(v) for v in s]
    x       = [pt[1] for pt in pts]
    y       = [pt[2] for pt in pts]
    return x,y
end

@recipe function f(els::NTuple{<:Any,<:AbstractElement{ReferenceLine}};npts=10)
    label --> ""
    grid   --> false
    for el in els
        @series begin
            el
        end
    end
end

@recipe function f(ent::AbstractElement{ReferenceSquare};h=0.1)
    legend --> false
    grid   --> true
    seriesalpha --> 0.1
    # aspect_ratio --> :equal
    seriestype --> :wireframe
    xrange = 0:h:1
    yrange = 0:h:1
    pts    = [ent((x,y)) for x in xrange, y in yrange]
    x      =  [pt[1] for pt in pts]
    y      =  [pt[2] for pt in pts]
    z      =  [pt[3] for pt in pts]
    x,y,z
end

@recipe function f(el::ParametricElement)
    sz = 10
    D = el.preimage
    grid   --> false
    aspect_ratio --> :equal
    label --> ""
    if D isa HyperRectangle{1}
        s       = LinRange(0,1,sz)
        pts     = [el(v) for v in s]
        x       = [pt[1] for pt in pts]
        y       = [pt[2] for pt in pts]
        x,y
    elseif D isa HyperRectangle{2}
        seriestype := :surface
        xrange = LinRange(0,1,sz)
        yrange = LinRange(0,1,sz)
        pts    = [el((x,y)) for x in xrange, y in yrange]
        x      =  [pt[1] for pt in pts]
        y      =  [pt[2] for pt in pts]
        z      =  [pt[3] for pt in pts]
        x,y,z
    else
        notimplemented()
    end
end

@recipe function f(els::Vector{<:ParametricElement})
    for el in els
        @series begin
            el
        end
    end
end
