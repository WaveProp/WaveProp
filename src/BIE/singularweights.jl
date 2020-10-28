function singular_weights(op::IntegralOperator,method=:kress)
    X,Y = op.X, op.Y
    dim = geometric_dimension(Y)
    near_list = near_interaction_list(X,Y;dim)
    return near_list
end    

function near_interaction_list(X,Y;dim)
    dict = Dict{DataType,Vector{Vector{Int}}}()    
    for E in etypes(Y)
        geometric_dimension(E) == dim || continue
        ilist = etype_near_interaction_list(E,X,Y)
        push!(dict,E=>ilist)
    end
    return dict
end    

function etype_near_interaction_list(E,X,Y)
    ilist = Vector{Vector{Int}}()
    iter  = ElementIterator{E}(Y)
    for el in ElementIterator{E}(Y)
        inear = element_near_interaction_list(el,X)        
        push!(ilist,inear)
    end 
    ilist       
end    

function element_near_interaction_list(el,X,tol=1e-2)
    ilist = Int[]
    vtx = el.nodes
    for (i,x) in enumerate(qnodes(X))
        for y in vtx
            d = norm(x-y)    
            if d < tol
                push!(ilist,i)    
                break # one point is already close, so move onto the next x
            end
        end
    end            
    return ilist
end    


