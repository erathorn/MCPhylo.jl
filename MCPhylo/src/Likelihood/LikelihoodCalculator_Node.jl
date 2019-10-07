
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node}, pi_::Number, rates::Vector{Float64})::Float64
    n_c = size(tree_postorder[1].data)[2]
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
    @simd for ind in 1:n_c
        @inbounds res += exp(log(rdata[1,ind])+ _lpi_) + exp(log(rdata[2,ind])+ _pi_)
    end

    return res#sum(log.(rdata.*[pi_, 1.0-pi_]))
end # function


function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Nothing
    #@assert size(node.child)[1] == 2
    @assert size(rates)[1] == n_c
    left_daughter::Node = node.lchild
    right_daughter::Node = node.rchild
    linc::Float64 = left_daughter.inc_length
    rinc::Float64 = right_daughter.inc_length
    left_daughter_data::Array{Float64,2} = left_daughter.data
    right_daughter_data::Array{Float64,2} = right_daughter.data

    # use the inbounds decorator to enable SIMD
    # SIMD greatly improves speed!!!
    Base.Threads.@threads for ind=eachindex(rates)
        @inbounds r::Float64 = rates[ind]

        @fastmath ext::Float64 = exp(-linc*r)
        ext_::Float64 = 1.0-ext
        p_::Float64 = 1.0-pi_
        v_::Float64 = ext_*pi_
        w_::Float64 = ext_*p_
        @fastmath v1::Float64 = log(ext+v_)
        @fastmath v2::Float64 = log(ext+w_)


        @inbounds a::Float64 = exp(log(left_daughter_data[1,ind])+v1) + exp(log(left_daughter_data[2,ind])+log(v_))
        @inbounds b::Float64 = exp(log(left_daughter_data[1,ind])+log(w_)) + exp(log(left_daughter_data[2,ind])+v2)

        @fastmath ext = exp(-rinc*r)
        ext_ = 1.0-ext
        v_ = ext_*pi_
        w_ = ext_*p_
        @fastmath v1 = log(ext+v_)
        @fastmath v2 = log(ext+w_)


        @inbounds c::Float64 = exp(log(right_daughter_data[1,ind])+v1) + exp(log(right_daughter_data[2,ind])+log(v_))
        @inbounds d::Float64 = exp(log(right_daughter_data[1,ind])+log(w_)) + exp(log(right_daughter_data[2,ind])+v2)

        @inbounds node.data[1,ind] = exp(log(a)+log(c))
        @inbounds node.data[2,ind] = exp(log(b)+log(d))
    end # for
end # function

function GradiantLog(tree_preorder::Vector{Node}, pi_::Number, rates::Array{Float64,1})
    root::Node = tree_preorder[1]
    n_c::Int64 = size(root.data)[2]
    Up::Array{Float64,3} = ones(length(tree_preorder)+1, size(root.data)[1], n_c)
    Grad_ll::Array{Float64,1} = zeros(length(tree_preorder))
    gradient::Array{Float64,1} = zeros(n_c)
    for node in tree_preorder
        node_ind::Int64 = node.num
        if node.binary == "1"
            # this is the root
            @simd for i in 1:n_c
                @inbounds Up[node_ind,1,i] = pi_
                @inbounds Up[node_ind,2,i] = 1.0-pi_
            end # for
        else
            sister::Node = get_sister(node)
            mother::Node = get_mother(node)


            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], sister.data, n_c)
            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], Up[mother.num,:,:], n_c)

            Base.Threads.@threads for i in 1:n_c
                @inbounds my_mat::Array{Float64,2} = exponentiate_binary(pi_, node.inc_length, rates[i])

                @inbounds a::Float64 = node.data[1, i] * my_mat[1,1] + node.data[2,i] * my_mat[2,1]
                @inbounds b::Float64 = node.data[1, i] * my_mat[1,2] + node.data[2,i] * my_mat[2,2]

                @inbounds gradient[i] = Up[node_ind,1,i] * a + Up[node_ind,2,i]*b

                @inbounds Up[node_ind,1,i] = Up[node_ind,1, i] * my_mat[1,2] + Up[node_ind,2,i] * my_mat[2,2]
                @inbounds Up[node_ind,2,i] = Up[node_ind,1, i] * my_mat[1,1] + Up[node_ind,2,i] * my_mat[2,1]
            end

            d = sum(pointwise_mat(Up[node_ind,:,:],node.data, n_c), dims=1)
            @inbounds gradient ./= d[1,:]
            @inbounds Grad_ll[node_ind] = sum(gradient)

            if node.nchild == 0
                @inbounds scaler = sum(Up[node_ind,:,:], dims=1)
                @inbounds Up[node_ind,:,:] ./= scaler
            end # if
        end # if
    end # for
    return Grad_ll

end # function
