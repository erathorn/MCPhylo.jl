
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

function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Array{Float64}
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
    Up::Array{Float64,3} = zeros(length(tree_preorder)+1, size(root.data)[1], size(root.data)[2])
    Grad_ll::Array{Float64} = zeros(length(tree_preorder))
    for node in tree_preorder
        if node.binary == "1"
            # this is the root
            Up[1,1,:] = pi_
            Up[1,2,:] = 1.0-pi_
        else
            sister::Node = Tree_Module.get_sister(root, node)
            node_ind::Int = parse(Int, node.binary, base=2)
            sister_ind::Int = parse(Int, sister.binary, base=2)
            Up[node_ind,:,:].*=Up[sister_ind,:,:]
            Up[node_ind,:,:].*=Up[node_ind,:,:]

            my_mat::Array{Float64,2} = exponentiate_binary(pi_, node.inc_length, r)

            a::Array{Float64,1} = node.data[1,:].*my_mat[1,1] .+ node.data[2,:].*my_mat[2,1]
            b::Array{Float64,1} = node.data[1,:].*my_mat[1,2] .+ node.data[2,:].*my_mat[2,2]

            gradient::Array{Float64,1} = Up[node_ind,1,:].*a .+ Up[node_ind,2,:].b

            Up[node_ind,1,:] = Up[node_ind,1,:].*my_mat[1,2] + Up[node_ind,2,:].*my_mat[2,2]
            Up[node_ind,2,:] = Up[node_ind,1,:].*my_mat[1,1] + Up[node_ind,2,:].*my_mat[2,1]

            d = sum(Up[node_ind,:,:].*node.data, dims=1)
            gradient ./= d
            Grad_ll[node_ind] = sum(gradient)

            if node.nchild == 0
                scaler = sum(Up[node_ind,:,:], dims=1)
                Up[node_ind,:,:] ./= scaler
            end # if
        end # if
    end # for
    return Grad_ll

end # function
