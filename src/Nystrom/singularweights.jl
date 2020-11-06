function _singular_weights_dim(op::IntegralOperator)
    X,Y = op.X, op.Y
    dim = geometric_dimension(Y)
    xs  = _source_location(Y)
    γ₀Ψ, γ₁Ψ = _traces(op,xs)
    R   = _residue()
end    

function singular_weights(op::IntegralOperator,method=:kress)
    X,Y = op.X, op.Y
    dim = geometric_dimension(Y)
    near_list = near_interaction_list(X,Y;dim=dim)
    _singular_weights_kress(op::IntegralOperator,near_list)
end 

function _singular_weights_kress(op,near_list)
    weights = []    
    jj      = []    
    X,Y = op.X, op.Y
    for (E,ilist) in near_list
        el_iter = ElementIterator{E}(Y)
        @assert length(ilist) == length(el_iter)
        nel = 1
        for (el,idx_near) in zip(el_iter,ilist)
            w = _singular_weights_kress(op,el,idx_near)
            j_global = el2qnodes(Y,E)[nel]
            push!(weights,w)
            push!(jj,j_global)
            nel += 1
        end    
    end 
    weights,jj
end    

function _singular_weights_kress(op,el,idx_near)
    X,Y = op.X, op.Y    
    K = op.kernel    
    for (i,jloc) in idx_near
        vi   = qnodes(qrule)
        x    = qnodes(op.X)[i]
        vs   = vi[jloc]
        # get the weights on the parametric nodes vi for the kernel constructed,
        # where `K` can be singular or nearly singular at el(vs). Note that we
        # operate solely in parameter space. 
        w    = qsing(vi,vs,el) do v
            K(x,el(v))*measure(el,v)    
        end
    end        
end    
 


