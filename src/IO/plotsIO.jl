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