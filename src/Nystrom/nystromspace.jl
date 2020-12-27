abstract type ApproximationSpace end

Base.@kwdef struct Quadrature{N,T}
    nodes::Vector{SVector{N,T}} = SVector{N,T}[]
    weights::Vector{T} = T[]
    normals::Vector{SVector{N,T}} = SVector{N,T}[]
end    

struct NystromElement{E,Q,Np,T}
    element::E
    quadrature::Q
    vals::SVector{Np,T}
end 

element(e::NystromElement)    = e.element
quadrature(e::NystromElement) = e.quadrature

"""
    struct Nystrom{N,T} <: ApproximationSpace
    
An `ApproximationSpace` given by a mesh and quadrature rule for each element
type in the mesh. 

The degrees of freedom of a `NystromSpace` are associated with its quadrature
nodes. Elements of this approximation space are of type `NystromElement`, and
can be constructed on-the-fly by combining an element of `mesh` with the
quadrature rule found on `etype2qrule`.

Use `elt2dofs[E][:,n]` to recover the degrees of freedom
associated with the `n`-th element of type `E`. 
"""
Base.@kwdef struct NystromSpace{N,T} <: ApproximationSpace
    mesh::GenericMesh{N,T}
    domain::Domain    
    etype2qrule::Dict{DataType,AbstractQuadratureRule}
    quadrature::Quadrature{N,T} = Quadrature{N,T}()
    elt2dof::Dict{DataType,Matrix{Int}} = Dict{DataType,Matrix{Int}}()
end  

function NystromSpace(mesh::GenericMesh,Ω::Domain,etype2qrule,compute_normal::Bool=true)
    N,T        = ambient_dimension(mesh), eltype(mesh)
    # initialize empty fields
    nys_mesh   = NystromSpace{N,T}(;domain=Ω,mesh=mesh,etype2qrule=etype2qrule)
    # loop over entities
    for ent in entities(Ω) 
        dict    = Dict{DataType,Vector{Int}}()    
        submesh = view(mesh,ent)
        for E in etypes(submesh)
            haskey(etype2qrule,E) || error("no quadrature rule found for element of type $E")
            qrule = etype2qrule[E]
            _build_nystrom_space!(nys_mesh,submesh,E,qrule,compute_normal)
            dict[E] = collect(istart:iend)
        end  
        nys_mesh.ent2elt[ent] = dict # new entry 
    end   
    return nys_mesh 
end    

@noinline function _build_nystrom_space!(nys_space,submesh,E,qrule,compute_normal)
    elt2dof         = get!(nys_mesh.elt2dof,E,Matrix{Int}(undef,0,0))
    elt2dof         = vec(elt2dof)
    x̂,ŵ             = qrule() #nodes and weights on reference element
    for el in elements(submesh,E)
        # compute quadrature on element
        x = map(el,x̂)
        w = map(zip(x̂,ŵ)) do (x̂,ŵ)
            ŵ*measure(el,x̂)    
        end 
        # append quadrature information
        istart = length(qnodes(nys_space)) + 1
        append!(qweights(nys_space),w)
        append!(qnodes(nys_space),x)
        if compute_normal
            ν = map(x->normal(el,x),x̂)   
            append!(qnormals(nys_space),ν)
        end
        iend   = length(qnodes(nys_mesh))
        # add information to recover qnodes from index of el
        append!(elt2dof,collect(istart:iend))
    end
    elt2dof = reshape(elt2dof,length(x̂),:)
    nys_space.el2dof[E] = el2dof
    return nys_space
end

dof(sp::NystromSpace) = 1:length(quadrature(sp))
qnodes(sp::NystromSpace)   = qnodes(quadrature(sp))
qweights(sp::NystromSpace) = qweights(quadrature(sp))
qnormals(sp::NystromSpace) = qnormals(quadrature(sp))


struct NystromSubSpace
    parent::NystromSpace
    domain::Domain
end    

# Σ        = NystromSpace

# Γ12_mesh = view(NystromSpace,Γ12)
# Γ13_mesh = view(NystromSpace,Γ13)

# S12 = SingleLayerOperator(op,Γ12_mesh,Σ)
# S21 = SingleLayerOperator(op,Γ21_mesh,Σ)

# mat = [S12;S21;S23]

