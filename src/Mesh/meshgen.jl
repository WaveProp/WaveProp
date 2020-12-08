function meshgen(p::AbstractParametricBody;gridsize)
    N       = ambient_dimension(p)    
    mesh    = GenericMesh{N,Float64}()
    dim,tag = Geometry._global_add_entity(p) # add entity to global list to make it easily retriavable using only (dim,tag)
    ent     = ElementaryEntity(dim,tag)
    Ω       = Domain(ent) 
    mesh.ent2tags[ent] = Dict{DataType,Vector{Int}}() # # no mesh is generated for domain, only its boundary, empty entry
    for bnd in boundary(p)
        # add elements for each boundary segment    
        bnd_dim,bnd_tag = Geometry._global_add_entity(bnd)
        bnd_ent         = ElementaryEntity(bnd_dim,bnd_tag)
        els             = _meshgen(bnd;gridsize)
        E               = eltype(els)
        push!(etypes(mesh),E)
        vals            = get!(mesh.els,E,Vector{E}())
        istart          = length(vals) + 1
        append!(vals,els)
        iend            = length(vals)
        # add mapping from entity to elements
        push!(boundary(ent),bnd_ent)
        mesh.ent2tags[bnd_ent] = Dict(E=>collect(istart:iend)) # add key
    end 
    unique!(etypes(mesh))   
    return Ω,mesh
end    

function _meshgen(p::ParametricEntity;gridsize)
    # mesh the domain of p
    f    = Geometry.parametrization(p)
    d    = domain(p)
    grid = CartesianMesh(d;gridsize)
    els = [Geometry.ParametricElement(f,d) for d in elements(grid)]
    return els
end    

# function meshgen!(mesh::GenericMesh,p::ParametricEntity,h)
#     # create a new ElementaryEntity for p    
#     dim  = geometric_dimension(p)    
#     ent  = ElementaryEntity(dim)    
#     # mesh the domain of p
#     f    = parametrization(p)
#     d    = domain(p)
#     grid = CartesianMesh(d;h=h)
#     # prepare the dictionary to keep track of the elements forming the entity to
#     # be meshed
#     T    = eltype(p)
#     F    = typeof(f)
#     E    = ParametricElement{D,T,F}
#     ent2tags = mesh.ent2tags
#     ent2tags[ent] = E => Int[]
#     # get vector storing elements of type E
#     if !haskey(mesh.els,E)
#         mesh.els[E] = Vector{E}()
#     end    
#     els  = mesh.els[E]
#     tags = ent2tags[ent][E]
#     # loop over the reference mesh and push the elements
#     tag_start = length(els)+1
#     for ref in elements(grid)
#         el = E(f,ref) # E = ParametricElement{D,T,F}
#         push!(els,el)
#     end    
#     tag_end  = length(els)
#     new_tags = collect(tag_start:tag_end)
#     append!(tags,new_tags)
#     return mesh
# end    

# function read_psurf(surf)
#     N    = ambient_dimension(surf)
#     mesh = GenericMesh{N,Float64}()
#     dim  = geometric_dimension(surf)
#     tag  = Geometry._new_tag(dim)
#     ent  = ElementaryEntity(dim)
#     Ω    = Domain([ent]) # empty domain
#     for bnd in boundary(surf)
#         dim_bnd = dim-1    
#         bnd = ElementaryEntity(dim_bnd)
#         push!(ent.boundary,bnd)
#     end
#     return Ω
# end   

# function meshgen!(ent::AbstractEntity{N,M},h) where {N,M}
#     domain   = ent.domain
#     n        = ceil.(Int,(domain.widths)/h)
#     elements = ent.elements
#     resize!(elements,prod(n))
#     if M == 1
#         x = range(domain.origin[1],domain.origin[1]+domain.widths[1],length=n[1]+1)
#         for i in 1:n[1]
#             elements[i] = HyperRectangle(x[i],x[i+1]-x[i])
#         end
#     elseif M == 2
#         x = range(domain.origin[1],domain.origin[1]+domain.widths[1],length=n[1]+1)
#         y = range(domain.origin[2],domain.origin[2]+domain.widths[2],length=n[2]+1)
#         cc = 1
#         for i in 1:n[1]
#             for j in 1:n[2]
#                 elements[cc] = HyperRectangle(x[i],y[j],x[i+1]-x[i],y[j+1]-y[j])
#                 cc += 1
#             end
#         end
#     else
#         error("type parameter $M must be 1 or 2")
#     end
#     return ent
# end

