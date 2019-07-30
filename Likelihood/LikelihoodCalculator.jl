
# TODO: Look into Parallelizing it
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node}, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Float64

    for node in tree_postorder
        if node.nchild !== 0
            CondLikeInternal(node, pi_, rates, n_c)
        end # if
    end # for

    # sum the two rows
    rdata::Array{Float64,2}=last(tree_postorder).data
    res::Float64 = 0.0
    _pi_::Float64 = log(1.0-pi_)
    _lpi_::Float64 = log(pi_)
    @inbounds for ind in 1:n_c
        res +=(log(rdata[1,ind])+ _lpi_) + (log(rdata[2,ind])+ _pi_)

        #rdata[1, ind] *pi_
        #rdata[2, ind] *=_pi_
    end

    return res#sum(log.(rdata.*[pi_, 1.0-pi_]))
end # function

"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree::Array{Float64,2}, data::Array{Float64,3}, mypi_::Number, rates::Vector{Float64}, n_c::Int64)::Float64
    tree_postorder::Vector{Int64} = post_order(tree)
    leaves::Vector{Int64} = get_leaves(tree)
    for node in tree_postorder
        if !(node in leaves)
            CondLikeInternal(tree, node, data, mypi_, rates, n_c)
        end # if
    end # for
    root::Int64 = find_root(tree)
    # sum the two rows
    #rdata::Array{Float64,2}=data[root,:,:]
    res::Float64 = 0.0
    _pi_::Float64 = log(1.0-mypi_)
    _lpi_::Float64 = log(mypi_)
    @inbounds for ind in 1:n_c
        res +=(log(data[1,ind, root])+ _lpi_) + (log(data[2,ind, root])+ _pi_)

        #rdata[1, ind] *pi_
        #rdata[2, ind] *=_pi_
    end

    return res#sum(log.(rdata.*[pi_, 1.0-pi_]))
end # function

function CondLikeInternal(tree::Array{Float64}, node::Int64, data::Array{Float64}, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Nothing
    @assert size(rates)[1] == n_c
    children::Vector{Int64} = get_neighbours(tree[node,:])
    left_daughter::Int64 = children[1]
    right_daughter::Int64 = children[2]
    linc::Float64 = tree[node, left_daughter]
    rinc::Float64 = tree[node, right_daughter]
    #@inbounds left_daughter_data::Array{Float64,2} = data[1:2, 1:n_c, left_daughter]
    #@inbounds right_daughter_data::Array{Float64,2} = data[1:2, 1:n_c, right_daughter]

    # use the inbounds decorator to enable SIMD
    # SIMD greatly improves speed!!!
    @inbounds for ind in 1:n_c
        @inbounds r::Float64 = rates[ind]
        left_mat::Array{Float64,2} = exponentiate_binary(pi_, linc, r)
        right_mat::Array{Float64,2} = exponentiate_binary(pi_, rinc, r)

        @inbounds a::Float64 = data[1,ind, left_daughter]*left_mat[1,1] + data[2,ind, left_daughter]*left_mat[2,1]
        @inbounds b::Float64 = data[1,ind, left_daughter]*left_mat[1,2] + data[2,ind, left_daughter]*left_mat[2,2]
        @inbounds c::Float64 = data[1,ind, right_daughter]*right_mat[1,1] + data[2,ind, right_daughter]*right_mat[2,1]
        @inbounds d::Float64 = data[1,ind, right_daughter]*right_mat[1,2] + data[2,ind, right_daughter]*right_mat[2,2]

        @inbounds data[1,ind, node] = a*c
        @inbounds data[2,ind, node] = b*d

    end # for
end # function




function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Nothing
    @assert size(node.child)[1] == 2
    @assert size(rates)[1] == n_c
    left_daughter::Node = node.child[1]
    right_daughter::Node = node.child[2]
    linc::Float64 = left_daughter.inc_length
    rinc::Float64 = right_daughter.inc_length
    left_daughter_data::Array{Float64,2} = left_daughter.data
    right_daughter_data::Array{Float64,2} = right_daughter.data

    # use the inbounds decorator to enable SIMD
    # SIMD greatly improves speed!!!
    @inbounds for ind in 1:n_c
        r::Float64 = rates[ind]
        left_mat::Array{Float64,2} = exponentiate_binary(pi_, linc, r)
        right_mat::Array{Float64,2} = exponentiate_binary(pi_, rinc, r)

        a::Float64 = left_daughter_data[1,ind]*left_mat[1,1] + left_daughter_data[2,ind]*left_mat[2,1]
        b::Float64 = left_daughter_data[1,ind]*left_mat[1,2] + left_daughter_data[2,ind]*left_mat[2,2]
        c::Float64 = right_daughter_data[1,ind]*right_mat[1,1] + right_daughter_data[2,ind]*right_mat[2,1]
        d::Float64 = right_daughter_data[1,ind]*right_mat[1,2] + right_daughter_data[2,ind]*right_mat[2,2]

        node.data[1,ind] = a*c
        node.data[2,ind] = b*d

    end # for
end # function

function GradiantLog(tree_preorder::Vector{Node}, pi_::Number)
    root::Node = tree_preorder[1]
    n_c::Int64 = size(root.data)[2]
    Up::Array{Float64,3} = ones(length(tree_preorder)+1, size(root.data)[1], n_c)
    Grad_ll::Array{Float64} = zeros(length(tree_preorder))
    for node in tree_preorder
        if node.binary == "1"
            # this is the root
            @inbounds for i in 1:n_c
                Up[node.num,1,i] = pi_
                Up[node.num,2,i] = 1.0-pi_
            end # for
        else
            sister::Node = get_sister(root, node)
            mother::Node = get_mother(root, node)
            node_ind::Int64 = node.num

            Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], sister.data, n_c)
            Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], Up[mother.num,:,:], n_c)
            #Up[node_ind,:,:].*=sister.data
            #Up[node_ind,:,:].*=Up[parse(Int, mother.binary, base=2),:,:]

            my_mat::Array{Float64,2} = exponentiate_binary(pi_, node.inc_length, 1.0)

            #a::Array{Float64,1} = node.data[1,:].*my_mat[1,1] .+ node.data[2,:].*my_mat[2,1]
            #b::Array{Float64,1} = node.data[1,:].*my_mat[1,2] .+ node.data[2,:].*my_mat[2,2]
            a::Array{Float64,1} = my_dot(node.data, my_mat[:,1], n_c)
            b::Array{Float64,1} = my_dot(node.data, my_mat[:,2], n_c)

            #gradient::Array{Float64,1} = Up[node_ind,1,:].*a .+ Up[node_ind,2,:].*b
            gradient::Array{Float64,1} = pointwise_vec(Up[node_ind,1,:],a, n_c) .+ pointwise_vec(Up[node_ind,2,:],b,n_c)

            #Up[node_ind,1,:] = Up[node_ind,1,:].*my_mat[1,2] + Up[node_ind,2,:].*my_mat[2,2]
            #Up[node_ind,2,:] = Up[node_ind,1,:].*my_mat[1,1] + Up[node_ind,2,:].*my_mat[2,1]
            Up[node_ind,1,:] = my_dot(Up[node_ind,:,:], my_mat[:,2], n_c)
            Up[node_ind,2,:] = my_dot(Up[node_ind,:,:], my_mat[:,1], n_c)

            #d = sum(Up[node_ind,:,:].*node.data, dims=1)
            d = sum(pointwise_mat(Up[node_ind,:,:],node.data, n_c), dims=1)
            gradient ./= d[1,:]
            Grad_ll[node_ind] = sum(gradient)

            if node.nchild == 0
                scaler = sum(Up[node_ind,:,:], dims=1)
                Up[node_ind,:,:] ./= scaler
            end # if
        end # if
    end # for
    return Grad_ll

end # function


function GradiantLog(tree::Array{Float64,2}, data::Array{Float64,3}, pi_::Number)
    root::Int64 = find_root(tree)
    tree_preorder = pre_order(tree)
    n_c::Int64 = size(data)[2]
    leaves::Vector = get_leaves(tree)

    Up::Array{Float64,3} = ones(length(tree_preorder)+1, size(data)[1], n_c)
    Grad_ll::Array{Float64} = zeros(length(tree_preorder))
    for node in tree_preorder
        if node==root
            # this is the root
            @inbounds for i in 1:n_c
                Up[node,1,i] = pi_
                Up[node,2,i] = 1.0-pi_
            end # for
        else
            mother::Int64 = get_mother(tree, node)
            n::Vector = get_neighbours(tree[mother,:])
            sister::Int64 = 0
            for i in n
                if i != node
                    sister = i
                end
            end


            Up[node,:,:] = pointwise_mat(Up[node,:,:], data[:,:,sister], n_c)
            Up[node,:,:] = pointwise_mat(Up[node,:,:], Up[mother,:,:], n_c)
            #Up[node_ind,:,:].*=sister.data
            #Up[node_ind,:,:].*=Up[parse(Int, mother.binary, base=2),:,:]

            my_mat::Array{Float64,2} = exponentiate_binary(pi_, tree[mother, node], 1.0)

            #a::Array{Float64,1} = node.data[1,:].*my_mat[1,1] .+ node.data[2,:].*my_mat[2,1]
            #b::Array{Float64,1} = node.data[1,:].*my_mat[1,2] .+ node.data[2,:].*my_mat[2,2]
            a::Array{Float64,1} = my_dot(data[:,:,node], my_mat[:,1], n_c)
            b::Array{Float64,1} = my_dot(data[:,:,node], my_mat[:,2], n_c)

            #gradient::Array{Float64,1} = Up[node_ind,1,:].*a .+ Up[node_ind,2,:].*b
            gradient::Array{Float64,1} = pointwise_vec(Up[node,1,:],a, n_c) .+ pointwise_vec(Up[node,2,:],b,n_c)

            #Up[node_ind,1,:] = Up[node_ind,1,:].*my_mat[1,2] + Up[node_ind,2,:].*my_mat[2,2]
            #Up[node_ind,2,:] = Up[node_ind,1,:].*my_mat[1,1] + Up[node_ind,2,:].*my_mat[2,1]
            Up[node,1,:] = my_dot(Up[node,:,:], my_mat[:,2], n_c)
            Up[node,2,:] = my_dot(Up[node,:,:], my_mat[:,1], n_c)

            #d = sum(Up[node_ind,:,:].*node.data, dims=1)
            d = sum(pointwise_mat(Up[node,:,:],data[:,:,node], n_c), dims=1)
            gradient ./= d[1,:]
            Grad_ll[node] = sum(gradient)

            if node in leaves
                scaler = sum(Up[node,:,:], dims=1)
                Up[node,:,:] ./= scaler
            end # if
        end # if
    end # for
    return Grad_ll

end # function
