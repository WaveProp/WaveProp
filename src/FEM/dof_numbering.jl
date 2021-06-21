abstract type DOFNumbering end

"""
Creates and stores the DOF numbering associated to the approximation space
constructed using the basis `b` in the `Domain` `Ω`.

In practice, the vector `dofs` contains all unique *global* identifiers (type
`T`) of degrees of freedom present in the domain `Ω` in reference to the
global information contained in `mesh`. The length of this vector depends on
the domain `Ω` but not the values of its elements. The *local* numbering of
the DOFs are given by the indices in this vector `dofs`, so in a sense this
gives a local to global numbering mapping.

The vector `elt2dof` contains, for each element (on which the basis `b` is
defined) in the domain `Ω`, the indices in the vector `dofs` of all degrees of
freedom that are (partly) supported on the said element.

Some optimization is possible in some cases where we could avoid the
allocation of one or both attributes `dofs` and `elt2dof` altogether.

Warning: this is not robust if the underlying `mesh` is modified.
"""
struct LocalDOFNumbering{S,T} <: DOFNumbering
    mesh::AbstractMesh
    Ω::Domain
    b::PolynomialBasis
    dofs::Vector{S}
    elt2dof::Vector{T}
end

dofs(dofnb::LocalDOFNumbering) = dofnb.dofs
elt2dof(dofnb::LocalDOFNumbering) = dofnb.elt2dof

function _element_support_basis(mesh, b::LagrangeBasis{D,p}) where {D,p}
    elts = DataType[]
    for E in etypes(mesh)
        if E <: AbstractElement
            if typeof(domain(E)) == D
                push!(elts, E)
            end
        end
    end
    @assert length(elts) == 1
    return elts[1]
end

"""
One degree of freedom per mesh element in the domain `Ω`.
"""
function local_dof_numbering(mesh::GenericMesh{d,T}, Ω::Domain,
                             b::LagrangeBasis{D,0}) where {d,D,Np,T}
    E = _element_support_basis(mesh, b)
    dofs = Int64[]
    for ent in entities(Ω)
        append!(dofs, ent2tags(mesh)[ent][E])
    end
    elt2dof = collect(1:length(dofs))
    return LocalDOFNumbering(mesh, Ω, b, dofs, elt2dof)
end

"""
One degree of freedom per mesh node in the domain `Ω`.
"""
function local_dof_numbering(mesh::GenericMesh, Ω::Domain,
                             b::LagrangeBasis{D,1}) where {D}
    E = _element_support_basis(mesh, b)
    N = number_of_nodes(D)
    node_tags = Int64[]
    for ent in entities(Ω)
        ent_tags = ent2tags(mesh)[ent][E]
        node_tags_ent = unique!(sort!(elements(mesh)[E][:,ent_tags][:]))
        append!(node_tags, node_tags_ent)
    end
    dofs = unique!(sort!(node_tags))
    glob2loc = SparseVector(length(nodes(mesh)),dofs,collect(1:length(dofs)))
    elt2dof = SVector{N,Int64}[]
    for ent in entities(Ω)
        ent_tags = ent2tags(mesh)[ent][E]
        node_tags_ent = elements(mesh)[E][:,ent_tags]
        for ie in 1:size(node_tags_ent,2)
            push!(elt2dof, glob2loc[node_tags_ent[:,ie]])
        end
    end
    return LocalDOFNumbering(mesh, Ω, b, dofs, elt2dof)
end

function assembly(mesh, Ω, u, v; order=1,
                  f=(u,v)->(i,j,x̂,_,_)->u(x̂)[j]*v(x̂)[i])
    submesh = view(mesh, Ω)
    etype2qrule = Mesh._qrule_for_mesh(submesh, order)
    E = _element_support_basis(submesh, u)
    @assert E == _element_support_basis(submesh, v)
    qrule = etype2qrule[E]
    udofnb = local_dof_numbering(mesh, Ω, u)
    vdofnb = local_dof_numbering(mesh, Ω, v)
    n = length(dofs(udofnb))
    m = length(dofs(vdofnb))
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    M  = zeros(Float64,length(v),length(u))
    for (iel, el) in enumerate(elements(submesh, E))
        utags = elt2dof(udofnb)[iel]
        vtags = elt2dof(vdofnb)[iel]
        elementary_matrix!(M,el, u, v, qrule; f=f)
        for (iloc, iglob) in enumerate(vtags)
            for (jloc, jglob) in enumerate(utags)
                push!(Is, iglob)
                push!(Js, jglob)
                push!(Vs, M[iloc,jloc])
            end
        end
        fill!(M,0)
    end
    A = sparse(Is,Js,Vs,m,n)
    return A
end

function assembly(mesh, Ω, v; order=1,
                  f=(v)->(i,x̂,_,_)->v(x̂)[i])
    submesh = view(mesh, Ω)
    etype2qrule = Mesh._qrule_for_mesh(submesh, order)
    E = _element_support_basis(submesh, v)
    qrule = etype2qrule[E]
    vdofnb = local_dof_numbering(mesh, Ω, v)
    m = length(dofs(vdofnb))
    b = zeros(Float64, m)
    for (iel, el) in enumerate(elements(submesh, E))
        vtags = elt2dof(vdofnb)[iel]
        M = elementary_matrix(el, v, qrule; f=f)
        for (iloc, iglob) in enumerate(vtags)
            b[iglob] += M[iloc]
        end
    end
    return b
end
