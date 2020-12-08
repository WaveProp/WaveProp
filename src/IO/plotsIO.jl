# plot a vector of points
@recipe function f(pts::AbstractVector{<:Point{N}}) where {N}
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

@recipe function f(pts::AbstractMatrix{<:Point}) 
    vec(pts)
end

# plot a hyperrectangle
@recipe function f(rec::HyperRectangle{N}) where {N}
    seriestype := :path
    label := ""
    if N == 2
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1],pt2[1]
        y1, y2 = pt1[2],pt2[2]
        [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1]
    elseif N == 3
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1],pt2[1]
        y1, y2 = pt1[2],pt2[2]
        z1, z2 = pt1[3],pt2[3]
        @series [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z1,z1,z1,z1,z1]
        @series [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z2,z2,z2,z2,z2]
        @series [x1,x1],[y1,y1],[z1,z2]
        @series [x2,x2],[y1,y1],[z1,z2]
        @series [x2,x2],[y2,y2],[z1,z2]
        @series [x1,x1],[y2,y2],[z1,z2]
    end
end

# recipe for parametric line
@recipe function f(ent::ParametricEntity{ReferenceLine};h=0.01)
    par = ent.parametrization
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    s       = 0:h:1    
    pts     = [par(v) for v in s]
    x       = [pt[1] for pt in pts]
    y       = [pt[2] for pt in pts]
    x,y
end

@recipe function f(ent::ParametricElement;gridsize=0.01)
    par = ent.parametrization
    legend --> false
    grid   --> false
    aspect_ratio --> :equal
    d = domain(ent)
    if d === ReferenceLine()
        h       =  gridsize[1]
        a       = ent.domain.low_corner[1]
        b       = ent.domain.high_corner[1]
        s       =  a:h:b
        pts     = [par(v) for v in s]
        x       = [pt[1] for pt in pts]
        y       = [pt[2] for pt in pts]
        return x,y
    else
        notimplemented()    
    end
end

# recipe for paramatric surface
@recipe function f(ent::ParametricEntity{ReferenceSquare};h=0.1)
    legend --> false
    grid   --> false
    aspect_ratio --> :equal
    seriestype := :surface
    xrange = 0:h:1
    yrange = 0:h:1
    pts    = [ent.parametrization((x,y)) for x in xrange, y in yrange]
    x      =  [pt[1] for pt in pts]
    y      =  [pt[2] for pt in pts]
    z      =  [pt[3] for pt in pts]
    vec(x),vec(y),vec(z)
end

# recipe for parametric body
@recipe function f(bdy::AbstractParametricBody)
    aspect_ratio --> :equal    
    for patch in boundary(bdy)
        @series begin
            patch
        end
    end
end

# recipe for many parametric bodies
@recipe function f(bdies::Vector{<:AbstractParametricBody})
    aspect_ratio --> :equal    
    for bdy in bdies
        @series begin
            bdy
        end
    end
end