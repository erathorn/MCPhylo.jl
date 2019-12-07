
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Any, pi_::Number, rates::Vector{Float64}, data::Array, n_c::Int64, blv)#::Float64
    mm::Array{Float64,2} = [0.0 0.0;0.0 0.0]
    mm2::Array{Float64,2} = [0.0 0.0;0.0 0.0]
    for node in tree_postorder
        if node.nchild != 0
            CondLikeInternal(node, pi_, rates, n_c, data, blv, mm, mm2)
        end
    end # for

    # sum the two rows
    rnum::Int64 = last(tree_postorder).num
    res = sum(log.(sum(data[:, : , rnum] .* Array([pi_, 1.0-pi_]), dims=1)))

    return res
end # function


"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Any, pi_::Number, rates::Vector{Float64}, data::CuArray, n_c::Int64, blv::Vector{Float64})#::Float64
    mm::CuArray{Float64} = [0 0;0 0]
    mm2::CuArray{Float64} = [0 0;0 0]

    for node in tree_postorder
        if node.nchild != 0
            CondLikeInternal(node, pi_, rates, n_c, data, blv)
        end
    end # for

    # sum the two rows
    rnum::Int64 = last(tree_postorder).num

    res = sum(CUDAnative.log.(sum(data[:, : , rnum] .* CuArray([pi_, 1.0-pi_]), dims=1)))

    return res
end # function


"""
This function calculates the log sum at an internal node in a tree. This function
is designed to specifically handle binary data. The matrix exponentiation is pulled
out of an extra function to avoid array creation.

The loop is written such that the julia converter should be able to infer simd patterns.
"""
function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64, data::Array, blv, mm, mm2)#::Float64

    @inbounds l_num::Int64 = node.lchild.num
    @inbounds r_num::Int64 = node.rchild.num
    @inbounds linc::Float64 = blv[l_num]
    @inbounds rinc::Float64 = blv[r_num]
    @inbounds n_num::Int64 = node.num

    r::Float64 = 1.0
    v::Float64 = exp(-r*linc)
    m11::Float64 = pi_ - (pi_ - 1)*v
    m22::Float64 = pi_*(v - 1) + 1
    m21::Float64 = pi_ - pi_* v
    m12::Float64 = (pi_ - 1)*(v - 1)
    mm[1,1] = m11
    mm[2,1] = m12
    mm[1,2] = m21
    mm[2,2] = m22
    v1::Float64 = exp(-r*rinc)
    m11b::Float64 = pi_ - (pi_ - 1)*v1
    m22b::Float64 = pi_*(v1 - 1) + 1
    m21b::Float64 = pi_ - pi_* v1
    m12b::Float64 = (pi_ - 1)*(v1 - 1)
    mm2[1,1] = m11b
    mm2[2,1] = m12b
    mm2[1,2] = m21b
    mm2[2,2] = m22b

    Base.Threads.@threads for ind in 1:n_c #eachindex(rates)

        @inbounds data[:,ind,n_num] = (mm*data[:,ind,l_num]).*(mm2*data[:,ind,r_num])
    end # for

end # function


function kernel(arr, trans_arr1, trans_arr2, n_num, l_num, r_num)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= size(arr,2)
        @inbounds arr[1,i,n_num]= ((arr[1,i,l_num]*trans_arr1[1,1]) + (arr[2,i,l_num]*trans_arr1[2,1])) *
                                  ((arr[1,i,r_num]*trans_arr2[1,1]) + (arr[2,i,r_num]*trans_arr2[2,1]))
        @inbounds arr[2,i,n_num]= ((arr[1,i,l_num]*trans_arr1[1,2]) + (arr[2,i,l_num]*trans_arr1[2,2])) *
                                  ((arr[1,i,r_num]*trans_arr2[1,2]) + (arr[2,i,r_num]*trans_arr2[2,2]))
    end
    return
end

function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64, data::CuArray, blv::Vector)

    @inbounds l_num::Int64 = node.lchild.num
    @inbounds r_num::Int64 = node.rchild.num
    @inbounds linc::Float64 = blv[l_num]
    @inbounds rinc::Float64 = blv[r_num]
    @inbounds n_num::Int64 = node.num

    r::Float64 = 1.0
    v::Float64 = exp(-r*linc)
    m11::Float64 = pi_ - (pi_ - 1)*v
    m22::Float64 = pi_*(v - 1) + 1
    m21::Float64 = pi_ - pi_* v
    m12::Float64 = (pi_ - 1)*(v - 1)
    mm::CuArray{Float64} = [m11 m21;m12 m22]

    v1::Float64 = exp(-r*rinc)
    m11b::Float64 = pi_ - (pi_ - 1)*v1
    m22b::Float64 = pi_*(v1 - 1) + 1
    m21b::Float64 = pi_ - pi_* v1
    m12b::Float64 = (pi_ - 1)*(v1 - 1)
    mm2::CuArray{Float64} = [m11b m21b;m12b m22]



    @cuda threads=1024 blocks=1 kernel(data, mm, mm2, n_num, l_num, r_num)


end # function





"""
    function GradiantLog(tree_preorder::Vector{Node}, pi_::Number, rates::Array{Float64,1}, data::Array{Fl0at64,3}, n_c::Int64)::Array{Float64}

This functio calculates the gradient for the Probabilistic Path Sampler.

"""

function GradiantLog(tree_preorder::Vector{Node}, pi_::Number, rates::Array{Float64,1},
                     data, n_c::Int64)::Array{Any}

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
            #sisterdata ./= sum(sisterdata, dims=1)

            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], sisterdata, n_c)
            @inbounds Up[node_ind,:,:] = pointwise_mat(Up[node_ind,:,:], Up[mother.num,:,:], n_c)

            nodedata::Array{Float64, 2} = exp.(data[node.num, :, :])
            #nodedata ./= sum(nodedata, dims=1)

            Base.Threads.@threads for i in 1:n_c
                @inbounds my_mat::Array{Float64,2} = exponentiate_binary_grad(pi_, node.inc_length, rates[i])
                @inbounds my_mat2::Array{Float64,2} = exponentiate_binary(pi_, node.inc_length, rates[i])

                @inbounds a = nodedata[1, i] * my_mat[1,1] + nodedata[2,i] * my_mat[1,2]
                @inbounds b = nodedata[1, i] * my_mat[2,1] + nodedata[2,i] * my_mat[2,2]

                @inbounds gradient[i] = Up[node_ind,1,i] * a + Up[node_ind,2,i]*b

                @inbounds Up[node_ind,1,i] = Up[node_ind,1, i] * my_mat2[1,1] + Up[node_ind,2,i] * my_mat2[2,1]
                @inbounds Up[node_ind,2,i] = Up[node_ind,1, i] * my_mat2[1,2] + Up[node_ind,2,i] * my_mat2[2,2]
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

function Grad_intern(blen_v, node_ind_v, binary_v, pi_, rates, data, n_c)

    for node_ind in node_ind_v
        #node_ind::Int64 = node.num

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
end
