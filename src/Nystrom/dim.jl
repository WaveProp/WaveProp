function assemble_dim(iop::IntegralOperator;compress=Matrix,location=:onsurface)
    X    = target_surface(iop)
    Y    = source_surface(iop)
    pde  = iop.kernel.pde
    T    = kernel_type(iop)
    if T === SingleLayer()
        singlelayer_dim(pde,X,Y;compress,location)
    elseif T === DoubleLayer()
        doublelayer_dim(pde,X,Y;compress,location)
    elseif T === AdjointDoubleLayer()
        adjointdoublelayer_dim(pde,X,Y;compress,location)
    elseif T === HyperSingular()
        hypersingular_dim(pde,X,Y;compress,location)
    else
        notimplemented()
    end
end

function single_doublelayer_dim(pde,X,Y=X;compress=Matrix,location=:onsurface)
    msg = "unrecognized value for kw `location`: received $location.
           Valid options are `:onsurface`, `:inside` and `:outside`."
    σ = location === :onsurface ? -0.5 : location === :inside ? 0 : location === :outside ? -1 : error(msg)
    Sop  = SingleLayerOperator(pde,X,Y)
    Dop  = DoubleLayerOperator(pde,X,Y)
    # convert to a possibly more efficient format
    S = compress(Sop)
    D = compress(Dop)
    basis,γ₁_basis = _basis_dim(Sop)
    γ₀B,γ₁B,R      = _auxiliary_quantities_dim(Sop,S,D,basis,γ₁_basis,σ)
    # compute corrections
    δS = _singular_weights_dim(Sop,γ₀B,γ₁B,R)
    δD = _singular_weights_dim(Dop,γ₀B,γ₁B,R)
    # add corrections to the dense part
    axpy!(true,δS,S)
    axpy!(true,δD,D)
    return S,D
end
singlelayer_dim(args...;kwargs...) = single_doublelayer_dim(args...;kwargs...)[1]
doublelayer_dim(args...;kwargs...) = single_doublelayer_dim(args...;kwargs...)[2]

function adjointdoublelayer_hypersingular_dim(pde,X,Y=X;compress=Matrix,location=:onsurface)
    msg = "unrecognized value for kw `location`: received $location.
    Valid options are `:onsurface`, `:inside` and `:outside`."
    σ = location === :onsurface ? -0.5 : location === :inside ? 0 : location === :outside ? -1 : error(msg)
    Kop  = AdjointDoubleLayerOperator(pde,X,Y)
    Hop  = HyperSingularOperator(pde,X,Y)
    # convert to a possibly more efficient compress
    K = compress(Kop)
    H = compress(Hop)
    basis,γ₁_basis = _basis_dim(Kop)
    γ₀B,γ₁B,R      = _auxiliary_quantities_dim(Kop,K,H,basis,γ₁_basis,σ)
    # compute corrections
    δK = _singular_weights_dim(Kop,γ₀B,γ₁B,R)
    δH = _singular_weights_dim(Hop,γ₀B,γ₁B,R)
    # add corrections to the dense part
    axpy!(true,δK,K)
    axpy!(true,δH,H)
    return K,H
end
adjointdoublelayer_dim(args...;kwargs...)  = adjointdoublelayer_hypersingular_dim(args...;kwargs...)[1]
hypersingular_dim(args...;kwargs...)       = adjointdoublelayer_hypersingular_dim(args...;kwargs...)[2]


function singular_weights_dim(iop::IntegralOperator,compress=Matrix)
    X,Y,op = iop.X, iop.Y, iop.kernel.op
    σ = X == Y ? -0.5 : 0.0
    #
    basis,γ₁_basis = _basis_dim(iop)
    Op1, Op2       = _auxiliary_operators_dim(iop,compress)
    γ₀B,γ₁B,R      = _auxiliary_quantities_dim(iop,Op1,Op2,basis,γ₁_basis,σ)
    Sp = _singular_weights_dim(iop,γ₀B,γ₁B,R)
    return Sp # a sparse matrix
end

function _auxiliary_quantities_dim(iop,Op1,Op2,basis,γ₁_basis,σ)
    T           = eltype(iop)
    X,Y    = target_surface(iop), source_surface(iop)
    nbasis = length(basis)
    # compute matrix of basis evaluated on Y
    ynodes   = qnodes(Y)
    γ₀B      = Matrix{T}(undef,length(ynodes),nbasis)
    γ₁B      = Matrix{T}(undef,length(ynodes),nbasis)
    for k in 1:nbasis
        for i in 1:length(ynodes)
            γ₀B[i,k] = basis[k](ynodes[i])
            γ₁B[i,k] = γ₁_basis[k](ynodes[i])
        end
    end
    # integrate the basis over Y
    xnodes   = qnodes(X)
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
                    R[i,k] += σ*γ₁_basis[k](xnodes[i])
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
    op = pde(kernel(iop))
    xs = _source_gen(iop)
    basis     = [(source) -> SingleLayerKernel(op)(x,source) for x in xs]
    γ₁_basis  = [(source) -> transpose(DoubleLayerKernel(op)(x,source)) for x in xs]
    return basis,γ₁_basis
end

function _singular_weights_dim(iop::IntegralOperator,γ₀B,γ₁B,R)
    X,Y = target_surface(iop), source_surface(iop)
    T   = eltype(iop)
    num_basis = size(γ₀B,2)
    a,b = combined_field_coefficients(iop)
    # we now have the residue R. For the correction we need the coefficients.
    dict_near = near_interaction_list(qnodes(X),Y;atol=0)
    Is = Int[]
    Js = Int[]
    Vs = T[]
    for (E,list_near) in dict_near
        el2qnodes = elt2dof(Y,E)
        num_qnodes, num_els   = size(el2qnodes)
        M                     = Matrix{T}(undef,2*num_qnodes,num_basis)
        @assert length(list_near) == num_els
        for n in 1:num_els
            j_glob                = @view el2qnodes[:,n]
            M[1:num_qnodes,:]     = @view γ₀B[j_glob,:]
            M[num_qnodes+1:end,:] = @view γ₁B[j_glob,:]
            # distinguish scalar and vectorial case
            if T <: Number
                F                     = qr(M)
            elseif T <: SMatrix
                F                     = qr!(blockmatrix_to_matrix(M))
            else
                error("unknown element type T=$T")
            end
            for (i,_) in list_near[n]
                if T <: Number
                    tmp = ((R[i:i,:])/F.R)*adjoint(F.Q)
                elseif T <: SMatrix
                    tmp_scalar  = (blockmatrix_to_matrix(R[i:i,:])/F.R)*adjoint(F.Q)
                    tmp  = matrix_to_blockmatrix(tmp_scalar,T)
                else
                    error("unknown element type T=$T")
                end
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
    Y      =  source_surface(iop)
    nquad  = 0
    for (E,tags) in elt2dof(Y)
        nquad = max(nquad,size(tags,1))
    end
    nbasis = 3*nquad
    # construct source basis
    return _source_gen(iop,nbasis;kfactor)
end

function _source_gen(iop,nsources;kfactor)
    N      = ambient_dimension(iop)
    Y      = source_surface(iop)
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

function _sphere_sources_lebedev(;nsources, radius=10, center=SVector(0.,0.,0.))
    lpts = lebedev_points(nsources)
    Xs = SVector{3,Float64}[]
    for pt in lpts
        push!(Xs,radius*pt .+ center)
    end
    return Xs
end

function _circle_sources(;nsources, radius=10, center=SVector(0.,0.))
    par   = (s) -> center .+ radius .* SVector(cospi(2 * s[1]), sinpi(2 * s[1]))
    x,_   = Integration._trapezoidalP(nsources)
    Xs    = SVector{2,Float64}[]
    for pt in x
        push!(Xs,par(pt))
    end
    return Xs
end
