function single_double_layer(pde,X,Y=X;compress=Matrix,correction=:dim,σ=-0.5)
    Sop  = SingleLayerOperator(pde,X,Y)
    Dop  = DoubleLayerOperator(pde,X,Y)
    # convert to a possibly more efficient format
    S = compress(Sop)
    D = compress(Dop)
    if correction == :dim
        basis,γ₁_basis = _basis_dim(Sop)    
        γ₀B,γ₁B,R      = _auxiliary_quantities_dim(Sop,S,D,basis,γ₁_basis,σ)
        # compute corrections
        δS = _singular_weights_dim(Sop,γ₀B,γ₁B,R)
        δD = _singular_weights_dim(Dop,γ₀B,γ₁B,R)
        # add corrections to the dense part
        axpy!(true,δS,S)
        axpy!(true,δD,D)
        return S,D
    elseif correction == :nothing
        return S,D
    else
        error("unrecognized correction method")
    end
end
singlelayer(args...;kwargs...) = single_double_layer(args...;kwargs...)[1]
doublelayer(args...;kwargs...) = single_double_layer(args...;kwargs...)[2]

