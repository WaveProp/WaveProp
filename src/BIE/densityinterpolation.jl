struct GreensCorrection{T,S,U} 
    iop::S
    R::Matrix{T}
    L::U
    idxel_near::Vector{Int}
end

Base.size(c::GreensCorrection,args...) = size(c.iop,args...)
Base.eltype(c::GreensCorrection{T}) where {T} = T

function GreensCorrection(iop::IntegralOperator,compress=Matrix)
    X,Y,op = iop.X, iop.Y, iop.kernel.op
    T = eltype(iop)
    # construct integral operators required for correction
    if kernel_type(iop) isa Union{SingleLayer,DoubleLayer}
        Op1 = IntegralOperator{T}(SingleLayerKernel(op),X,Y) |> compress
        Op2 = IntegralOperator{T}(DoubleLayerKernel(op),X,Y) |> compress
    elseif kernel_type(iop) isa Union{AdjointDoubleLayer,HyperSingular}
        Op1 = IntegralOperator{T}(AdjointDoubleLayerKernel(op),X,Y) |> compress
        Op2 = IntegralOperator{T}(HyperSingularKernel(op),X,Y) |> compress
    end
    GreensCorrection(iop,Op1,Op2)
end

function GreensCorrection(iop::IntegralOperator,Op1,Op2)
    xs = _source_gen(iop)
    GreensCorrection(iop,Op1,Op2,xs)
end

function GreensCorrection(iop::IntegralOperator,Op1,Op2,xs::Vector{<:Point})
    # construct greens "basis" from source locations xs
    op        = iop.kernel.op
    basis     = [y->SingleLayerKernel(op)(x,y) for x in xs]
    γ₁_basis  = [(y,ny)->transpose(DoubleLayerKernel(op)(x,y,ny)) for x in xs]
    #FIXME: the value of σ depends on whether the observation point in `X` lies
    #inside, outside, or on the integration surface. The current way of doing
    #is hacky: if the surfaces are the same we say the points are on, otherwise
    #we assume they are outside.
    σ = iop.X === iop.Y ? -0.5 : 0.0
    GreensCorrection(iop,Op1,Op2,basis,γ₁_basis,σ)
end

function GreensCorrection(iop::IntegralOperator,Op1,Op2,basis,γ₁_basis,σ)
    T           = eltype(iop)    
    kernel,X,Y  = iop.kernel, iop.X, iop.Y
    a,b         = combined_field_coefficients(kernel)
    op          = kernel.op
    m,n         = length(X),length(Y)
    nbasis      = length(basis)
    # compute matrix of basis evaluated on Y
    ynodes   = qnodes(Y)
    ynormals = qnormals(Y)
    γ₀B     = Matrix{T}(undef,length(ynodes),nbasis)
    γ₁B     = Matrix{T}(undef,length(ynodes),nbasis)
    for k in 1:nbasis
        for i in 1:length(ynodes)
            γ₀B[i,k] = basis[k](ynodes[i])
            γ₁B[i,k] = γ₁_basis[k](ynodes[i],ynormals[i])
        end
    end
    # integrate the basis over Y
    xnodes   = qnodes(X)
    xnormals = qnormals(X)
    R        = Op1*γ₁B - Op2*γ₀B   
    # analytic correction for on-surface evaluation of Greens identity
    if kernel_type(iop) isa Union{SingleLayer,DoubleLayer}
        if σ !== 0
            for k in 1:nbasis
                for i in 1:length(xnodes)
                    R[i,k] += σ*basis[k](xnodes[i])
                end
            end
        end
    elseif kernel_type(iop) isa Union{AdjointDoubleLayer,HyperSingular}
        if σ !== 0
            for k in 1:nbasis
                for i in 1:length(xnodes)
                    R[i,k] += σ*γ₁_basis[k](xnodes[i],xnormals[i])
                end
            end
        end
    end
    # we now have the residue R. For the correction we need the coefficients.
    dict_near = near_interaction_list(X,Y,dim=ambient_dimension(Y)-1,atol=1e-16)
    Is = Int[]
    Js = Int[]
    Vs = T[]
    for (E,list_near) in dict_near
        el2qnodes = Y.el2qnodes[E]
        num_qnodes, num_els   = size(el2qnodes)
        M                     = Matrix{T}(undef,2*num_qnodes,nbasis)
        @assert length(list_near) == num_els
        for n in 1:num_els
            j_glob             = el2qnodes[:,n]
            M[1:num_qnodes,:]     = γ₀B[j_glob,:]
            M[num_qnodes+1:end,:] = γ₁B[j_glob,:]
            F                  = qr(M)
            for (i,_) in list_near[n]
                tmp  = (R[i:i,:]/F.R)*adjoint(F.Q)
                w    = axpby!(a,view(tmp,1:num_qnodes),b,view(tmp,(num_qnodes+1):(2*num_qnodes)))
                append!(Is,fill(i,num_qnodes))
                append!(Js,j_glob)
                append!(Vs,w)
            end
        end    
    end        
    Sp = sparse(Is,Js,Vs,size(iop)...)
    return Sp
end

# FIXME:this function needs to be cleaned up and optimized for perf
function precompute_weights_qr(c::GreensCorrection)
    T = eltype(c)    
    w = [Vector{T}() for _ in 1:size(c,1)]
    iop  = c.iop
    a,b  = combined_field_coefficients(iop.kernel)
    X,Y  = iop.X, iop.Y
    QRTYPE = Base.promote_op(qr,Matrix{T})
    F      = Vector{QRTYPE}(undef,length(c.L))
    for i in 1:size(c,1)
        idx_el    = c.idxel_near[i]
        idx_el < 0 && continue
        idx_nodes = iop.Y.el2qnodes[idx_el]
        num_qnodes = length(idx_nodes)
        if !isassigned(F,idx_el)
            F[idx_el]  = qr(c.L[idx_el])
        end
        tmp     = (c.R[i:i,:]/F[idx_el].R)*adjoint(F[idx_el].Q)
        w[i]     = axpby!(a,view(tmp,1:num_qnodes),b,view(tmp,(num_qnodes+1):(2*num_qnodes)))
    end
    return w
end

function Base.Matrix(c::GreensCorrection)
    iop = c.iop
    a,b   = combined_field_coefficients(iop.kernel)
    M     = zeros(eltype(c),size(c))
    w     = precompute_weights_qr(c)
    for i in 1:size(c,1)
        idx_el    = c.idxel_near[i]
        idx_el < 0 && continue
        idx_nodes = iop.Y.el2qnodes[idx_el]
        M[i,idx_nodes] = w[i]
    end
    return M
end

function SparseArrays.sparse(c::GreensCorrection)
    m,n = size(c)
    iop = c.iop
    T   = eltype(c)
    a,b = combined_field_coefficients(iop.kernel)
    I = Int[]
    J = Int[]
    V = T[]
    w = precompute_weights_qr(c)
    for i in 1:size(c,1)
        idx_el    = c.idxel_near[i]
        idx_el < 0 && continue
        idx_nodes = getelements(iop.Y)[idx_el]
        append!(I,fill(i,length(idx_nodes)))
        append!(J,idx_nodes)
        append!(V,w[i])
    end
    return sparse(I,J,V,m,n)
end


function _source_gen(iop::IntegralOperator,kfactor=5)
    Y      =  iop.Y
    nquad  = mapreduce(x->length(x),max,Y.el2qnodes) # maximum number of quadrature nodes per element
    nbasis = 3*nquad
    # construct source basis
    return _source_gen(iop,nbasis,kfactor=kfactor)
end

function _source_gen(iop,nsources;kfactor)
    N      = ambient_dimension(iop)
    Y      = iop.Y
    pts    = qnodes(Y)
    # create a bounding box
    bbox   = bounding_box(pts)
    xc     = center(bbox)
    d      = diameter(bbox)
    if N == 2
        xs = _circle_sources(nsources=nsources,center=xc,radius=kfactor*d/2)
    elseif N == 3
        xs = _sphere_sources_lebedev(nsources=nsources,center=xc,radius=kfactor*d/2)
    else
        error("dimension must be 2 or 3. Got $N")
    end
    return xs
end

function _sphere_sources_lebedev(;nsources, radius=10, center=Point(0.,0.,0.))
    lpts = lebedev_points(nsources)
    Xs = Point{3,Float64}[]
    for pt in lpts
        push!(Xs,radius*pt .+ center)
    end
    return Xs
end