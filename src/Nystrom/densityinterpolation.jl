function singular_weights_dim(iop::IntegralOperator,compress=Matrix)
    X,Y,op = iop.X, iop.Y, iop.kernel.op        
    σ = X == Y ? -0.5 : 0.0
    # 
    basis,γ₁_basis = _basis_dim(iop)
    Op1, Op2       = _auxiliary_operators_dim(iop,compress)    
    γ₀B,γ₁B,R      = _auxiliary_quantities_dim(iop,Op1,Op2,basis,γ₁_basis,σ)
    _singular_weights_dim(iop,γ₀B,γ₁B,R)
end

function _auxiliary_quantities_dim(iop,Op1,Op2,basis,γ₁_basis,σ) 
    T           = eltype(iop)    
    kernel,X,Y  = iop.kernel, iop.X, iop.Y
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
    return γ₀B, γ₁B, R
end    

function _auxiliary_operators_dim(iop,compress)
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
    return Op1,Op2
end 

function _basis_dim(iop)
    op = iop.kernel.op
    xs = _source_gen(iop)
    basis     = [y->SingleLayerKernel(op)(x,y) for x in xs]
    γ₁_basis  = [(y,ny)->transpose(DoubleLayerKernel(op)(x,y,ny)) for x in xs]
    return basis,γ₁_basis
end    

function _singular_weights_dim(iop::IntegralOperator,γ₀B,γ₁B,R)
    X,Y = iop.X, iop.Y    
    T = eltype(iop)
    num_basis = size(γ₀B,2)
    a,b = combined_field_coefficients(iop)
    # we now have the residue R. For the correction we need the coefficients.
    dict_near = near_interaction_list(X,Y;dim=ambient_dimension(Y)-1,atol=1e-16)
    Is = Int[]
    Js = Int[]
    Vs = T[]
    for (E,list_near) in dict_near
        el2qnodes = Y.el2qnodes[E]
        num_qnodes, num_els   = size(el2qnodes)
        M                     = Matrix{T}(undef,2*num_qnodes,num_basis)
        @assert length(list_near) == num_els
        for n in 1:num_els
            j_glob                = el2qnodes[:,n]
            M[1:num_qnodes,:]     = γ₀B[j_glob,:]
            M[num_qnodes+1:end,:] = γ₁B[j_glob,:]
            F                     = qr(M)
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

function _source_gen(iop::IntegralOperator,kfactor=5)
    Y      =  iop.Y
    nquad  = 0
    for (E,tags) in el2qnodes(Y)        
        nquad = max(nquad,size(tags,1))
    end  
    nbasis = 3*nquad
    # construct source basis
    return _source_gen(iop,nbasis;kfactor)
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
        xs = _circle_sources(;nsources,center=xc,radius=kfactor*d/2)
    elseif N == 3
        xs = _sphere_sources_lebedev(;nsources,center=xc,radius=kfactor*d/2)
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

function _circle_sources(;nsources, radius=10, center=Point(0.,0.))
    geo   = Circle(center=center,radius=radius)
    par   = boundary(geo)[1]
    x,_   = Integration._trapezoidalP(nsources)
    Xs    = Point{2,Float64}[]
    for pt in x
        push!(Xs,par(pt))
    end
    return Xs
end
