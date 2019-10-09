
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node}, pi_::Number, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64


    for node in tree_postorder
        if node.nchild != 0
            CondLikeInternal(node, pi_, rates, n_c, data)
        end # if
    end # for

    # sum the two rows
    #rdata::Array{Float64,2}=last(tree_postorder).data
    rnum = last(tree_postorder).num
    res::Float64 = 0.0
    _lpi_::Float64 = log(1.0-pi_)
    _pi_::Float64 = log(pi_)

    @simd for ind in 1:n_c
        @inbounds res += log(exp(data[rnum,1,ind]+ _lpi_) + exp(data[rnum,2,ind]+ _pi_))
    end # for
    if res > 0
        blv = get_branchlength_vector(tree_postorder)
        println("res ", res, " blv ", blv, " pi ", pi_, " av ", unique(rates))
        @assert res < 0
    end
    return res
end # function


function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64, data::Array{Float64})::Nothing

    @assert size(rates)[1] == n_c
    left_daughter::Node = node.lchild
    right_daughter::Node = node.rchild
    linc::Float64 = left_daughter.inc_length
    rinc::Float64 = right_daughter.inc_length
    l_num::Int64 = left_daughter.num
    r_num::Int64 = right_daughter.num


    n_num::Int64 = node.num

    Base.Threads.@threads for ind=eachindex(rates)
        @inbounds r::Float64 = rates[ind]
        my_mat = log.(exponentiate_binary(pi_, linc, r))

        @inbounds a::Float64 = log(exp(data[l_num, 1,ind]+my_mat[1,1]) + exp(data[l_num, 2,ind]+my_mat[1,2]))
        @inbounds b::Float64 = log(exp(data[l_num, 1,ind]+my_mat[2,1]) + exp(data[l_num, 2,ind]+my_mat[2,2]))

        my_mat2 = log.(exponentiate_binary(pi_, rinc, r))

        @inbounds c::Float64 = log(exp(data[r_num, 1,ind]+my_mat2[1,1]) + exp(data[r_num,2,ind]+my_mat2[1,2]))
        @inbounds d::Float64 = log(exp(data[r_num, 1,ind]+my_mat2[2,1]) + exp(data[r_num,2,ind]+my_mat2[2,2]))
        if (a+c) > 0 || (b+d) > 0
            println("here ", a+c," ", b+d)
            println(a, ", ", c, ", ", b, ", ", d)
            println(data[l_num, 1,ind], " ", data[l_num, 2,ind])
            println(data[r_num, 1,ind], " ", data[r_num, 2,ind])
            println(my_mat)
            println(exp.(my_mat))
            println(my_mat2)
            println(exp.(my_mat2))
            throw("I am Here")
        end
        @inbounds data[n_num,1,ind] = a+c
        @inbounds data[n_num,2,ind] = b+d

    end # for
    if any(0 .< data[n_num, :, :])
        println(n_num)
        inds = findall(x-> x > 0 , data[n_num, :, :])
        for ind in inds
            println(data[n_num, :, :][ind])
        end
        @assert all(0 .>= data[n_num, :, :])
    end
end # function

function GradiantLog(tree_preorder::Vector{Node}, pi_::Number, rates::Array{Float64,1}, data, n_c)
    root::Node = tree_preorder[1]

    Up::Array{Float64,3} = ones(length(tree_preorder)+1, size(root.data)[1], n_c)
    Grad_ll::Array{Float64,1} = zeros(length(tree_preorder))
    gradient::Vector{Float64} = zeros(n_c)
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


            sisterdata::Array{Float64} = exp.(data[sister.num, :, :])
            sisterdata ./= sum(sisterdata, dims=1)
            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], sisterdata, n_c)


            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], Up[mother.num,:,:], n_c)

            nodedata::Array{Float64, 2} = exp.(data[node.num, :, :])
            nodedata ./= sum(nodedata, dims=1)

            Base.Threads.@threads for i in 1:n_c
                @inbounds my_mat::Array{Float64,2} = exponentiate_binary_grad(pi_, node.inc_length, rates[i])
                @inbounds my_mat2::Array{Float64,2} = exponentiate_binary(pi_, node.inc_length, rates[i])

                @inbounds a::Float64 = nodedata[1, i] * my_mat[1,1] + nodedata[2,i] * my_mat[1,2]
                @inbounds b::Float64 = nodedata[1, i] * my_mat[2,1] + nodedata[2,i] * my_mat[2,2]

                @inbounds gradient[i] = Up[node_ind,1,i] * a + Up[node_ind,2,i]*b

                @inbounds Up[node_ind,1,i] = Up[node_ind,1, i] * my_mat2[2,1] + Up[node_ind,2,i] * my_mat2[2,2]
                @inbounds Up[node_ind,2,i] = Up[node_ind,1, i] * my_mat2[1,1] + Up[node_ind,2,i] * my_mat2[1,2]
            end

            d = sum(pointwise_mat(Up[node_ind,:,:],nodedata, n_c), dims=1)

            @inbounds gradient ./= d[1,:]

            @inbounds Grad_ll[node_ind] = sum(gradient)

            if node.nchild != 0
                @inbounds scaler = sum(Up[node_ind,:,:], dims=1)
                @inbounds Up[node_ind,:,:] ./= scaler
            end # if
        end # if
    end # for

    return Grad_ll

end # function
