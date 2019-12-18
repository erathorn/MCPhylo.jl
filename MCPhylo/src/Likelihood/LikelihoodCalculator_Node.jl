
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node_ncu}, pi_::Number, rates::Vector{Float64}, data::Array, n_c::Int64, blv::Array{Float64,1})#::Float64
    #mm::Array{ForwardDiff.Dual} = [0.0 0.0;0.0 0.0]
    #mm2::Array{ForwardDiff.Dual} = [0.0 0.0;0.0 0.0]
    res = 0.0
    scaler = zeros(1, size(last(tree_postorder).data, 2))
    for node in tree_postorder
        if node.nchild != 0
            ld::Array{Float64, 2} = node.lchild.data
            rd::Array{Float64, 2} = node.rchild.data
            @inbounds l_num::Int64 = node.lchild.num
            @inbounds r_num::Int64 = node.rchild.num

            @inbounds linc::Float64 = blv[l_num]
            @inbounds rinc::Float64 = blv[r_num]
            out = CondLikeInternal(pi_, ld, rd, linc, rinc)
            m_s = maximum(out, dims = 1)
            out ./= m_s

            scaler .+= log.(m_s)
            node.data = out#CondLikeInternal(pi_, ld, rd, linc, rinc)
        end
    end # for

    # sum the two rows
    rnum::Array{Float64, 2} = last(tree_postorder).data#::Array{Number, 2}#.num

    res += sum(scaler.+log.(sum(rnum .* Array([pi_, 1.0-pi_]), dims=1)))#::Number

    return res
end # function


"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Array{Node_cu, 1}, pi_::Number, rates::Vector{Float64}, data::CuArray, n_c::Int64, blv::Array{Float64,1})#::Float64
    #mm::CuArray{Float64} = [0 0;0 0]
    #mm2::CuArray{Float64} = [0 0;0 0]
    p = convert(Float64, pi_)
    for node in tree_postorder
        if node.nchild != 0

            ld::CuArray{Float64, 2, Nothing} = node.lchild.data
            rd::CuArray{Float64, 2, Nothing} = node.rchild.data
            @inbounds l_num::Int64 = node.lchild.num
            @inbounds r_num::Int64 = node.rchild.num
            @inbounds linc::Float64 = blv[l_num]
            @inbounds rinc::Float64 = blv[r_num]
            v = exp(-linc)
            v1 = exp(-rinc)

            @cuda threads=128 blocks=16 kernel(node.data, p, v, v1, ld, rd)

        end
    end # for

    # sum the two rows
    rnum = last(tree_postorder).data::CuArray
    pia = CuArray([p, 1.0-p])
    res = CUDAnative.sum(CUDAnative.log.(CUDAnative.sum(rnum .* pia, dims=1)))::Float64

    return res
end # function




"""
This function calculates the log sum at an internal node in a tree. This function
is designed to specifically handle binary data. The matrix exponentiation is pulled
out of an extra function to avoid array creation.

The loop is written such that the julia converter should be able to infer simd patterns.
"""
function CondLikeInternal(pi_::Number, l_data::Array{Float64}, r_data::Array{Float64}, linc::Float64, rinc::Float64)#::Float64


        r::Float64 = 1.0
        mu =  1.0 / (2.0 * pi_ * (1-pi_))
        v = exp(-r*linc*mu)
        m11::Float64 = pi_ - (pi_ - 1)*v
        m22::Float64 = pi_*(v - 1) + 1
        m21::Float64 = pi_ - pi_* v
        m12::Float64 = (pi_ - 1)*(v - 1)
        mm = [m11 m12; m21 m22]

        v1 = exp(-r*rinc*mu)
        m11b::Float64 = pi_ - (pi_ - 1)*v1
        m22b::Float64 = pi_*(v1 - 1) + 1
        m21b::Float64 = pi_ - pi_* v1
        m12b::Float64 = (pi_ - 1)*(v1 - 1)
        mm2 = [m11b m12b; m21b m22b]

        #out = (mm*l_data) .* (mm2*r_data)
        out = (mm*l_data) .* (mm2*r_data)
        #out ./= sum(out, dims=1)
        return out #(mm*l_data) .* (mm2*r_data)  #[i.value for i in out]#(mm*l_data) .* (mm2*r_data)


end # function


function kernel(arr, pi_, v::Float64, v1::Float64, ldata, rdata)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= size(arr,2)
        m11 = pi_ - (pi_ - 1)*v
        m22 = pi_*(v - 1) + 1
        m21 = pi_ - pi_* v
        m12 = (pi_ - 1)*(v - 1)
        m11b = pi_ - (pi_ - 1)*v1
        m22b = pi_*(v1 - 1) + 1
        m21b = pi_ - pi_* v1
        m12b = (pi_ - 1)*(v1 - 1)
        @inbounds arr[1,i]= ((ldata[1,i]*m11) + (ldata[2,i]*m12)) * ((rdata[1,i]*m11b) + (rdata[2,i]*m12b))
        @inbounds arr[2,i]= ((ldata[1,i]*m21) + (ldata[2,i]*m22)) * ((rdata[1,i]*m21b) + (rdata[2,i]*m22b))
        #
    end
    return
end


function CondLikeInternal_c(pi_::Number, l_data::CuArray, r_data::CuArray, linc::Float64, rinc::Float64, arr::CuArray{Float64, 2})

    #r::Float64 = 1.0
    #v = exp(-r*linc)
    #m11::Float64 = pi_ - (pi_ - 1)*v
    #m22::Float64 = pi_*(v - 1) + 1
    #m21::Float64 = pi_ - pi_* v
    #m12::Float64 = (pi_ - 1)*(v - 1)
    #mm = CuArray([m11 m12; m21 m22])

    #v1 = exp(-r*rinc)
    #m11b::Float64 = pi_ - (pi_ - 1)*v1
    #m22b::Float64 = pi_*(v1 - 1) + 1
    #m21b::Float64 = pi_ - pi_* v1
    #m12b::Float64 = (pi_ - 1)*(v1 - 1)
    #mm2 = CuArray([m11b m12b; m21b m22b])

    @cuda threads = 1024 blocks = 32 kernel(arr, pi_, linc, rinc, l_data, r_data)
    #return (mm*l_data) .* (mm2*r_data) #[i.value for i in out]#(mm*l_data) .* (mm2*r_data)

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
