function singular_weights(op::IntegralOperator,method=:kress)
    X,Y = op.X, op.Y
    dim = geometric_dimension(Y)
    near_list = near_interaction_list(X,Y;dim)
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

function near_interaction_list(X,Y;dim)
    dict = Dict{DataType,Vector{Vector{Tuple{Int,Int}}}}()    
    for E in etypes(Y)
        geometric_dimension(E) == dim || continue
        ilist = etype_near_interaction_list(X,Y,E)
        push!(dict,E=>ilist)
    end
    return dict
end    

function etype_near_interaction_list(X,Y,E)
    ilist = Vector{Vector{Tuple{Int,Int}}}()
    e2n   = el2qnodes(Y,E)
    npts,nel = size(e2n)
    for n in 1:nel
        ynodes = qnodes(Y)[e2n[:,n]]
        inear  = _near_interaction_list(X,ynodes)
        push!(ilist,inear)
    end 
    ilist       
end    

function _near_interaction_list(X,ynodes,tol=1e-2)
    ilist    = Tuple{Int,Int}[]
    for (i,x) in enumerate(qnodes(X))            
        d,j = findmin([norm(x-y) for y in ynodes])        
        if d < tol
            push!(ilist,(i,j))    
        end
    end            
    return ilist
end    


