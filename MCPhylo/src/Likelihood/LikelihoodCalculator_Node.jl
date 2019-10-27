
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Any, pi_::Number, rates::Vector{Float64}, data::Array{Float64}, n_c::Int64, blv::Vector{Float64})#::Float64

    res = 0.0
    for node in tree_postorder
        if node.nchild != 0
            res += CondLikeInternal(node, pi_, rates, n_c, data, blv)
        end # if
    end # for

    # sum the two rows
    rnum = last(tree_postorder).num

    _lpi_::Float64 = log(1.0-pi_)
    _pi_::Float64 = log(pi_)

    @simd for ind in 1:n_c
        @inbounds res += log(exp(data[rnum,1,ind]+ _pi_) + exp(data[rnum,2,ind]+ _lpi_))
    end # for

    return res
end # function


"""
This function calculates the log sum at an internal node in a tree. This function
is designed to specifically handle binary data. The matrix exponentiation is pulled
out of an extra function to avoid array creation.

The loop is written such that the julia converter should be able to infer simd patterns.
"""
function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64, data::Array{Float64}, blv::Vector{Float64})::Float64

    @assert size(rates)[1] == n_c
    left_daughter::Node = node.lchild
    right_daughter::Node = node.rchild
    #linc::Float64 = left_daughter.inc_length
    #rinc::Float64 = right_daughter.inc_length
    l_num::Int64 = left_daughter.num
    r_num::Int64 = right_daughter.num
    linc = blv[l_num]
    rinc = blv[r_num]
    res = 0.0

    n_num::Int64 = node.num

    Base.Threads.@threads for ind=eachindex(rates)
        @inbounds r::Float64 = rates[ind]


        #ext = exp(-linc*r)
        #ext_ = 1.0-ext
        #pj_::Float64 = 1.0-pi_
        #m_v = ext+ext_*pi_
        #m_w = ext+ext_*pj_
        m = expo1(pi_, linc, r)
        if any(m .< 0)
            println(m)
            println(pi_, " ", linc, " ", r)
        end

        @inbounds a = log(exp(data[l_num, 1,ind]+log(m[1,1])) + exp(data[l_num, 2,ind]+log(m[2,1])))
        @inbounds b = log(exp(data[l_num, 1,ind]+log(m[1,2])) + exp(data[l_num, 2,ind]+log(m[2,2])))
        #my_mat2 = log.(exponentiate_binary(pi_, rinc, r))
        m = expo1(pi_, rinc, r)
        if any(m .< 0)
            println(m)
            println(pi_, " ", linc, " ", r)
        end
        #ext = exp(-rinc*r)
        #ext_ = 1.0-ext
        #v_ = ext_*pi_
        #w_ = ext_*pj_
        #m_v = ext+ext_*pi_
        #m_w = ext+ext_*pj_

        @inbounds c = log(exp(data[r_num, 1,ind]+log(m[1,1])) + exp(data[r_num,2,ind]+log(m[2,1])))
        @inbounds d = log(exp(data[r_num, 1,ind]+log(m[1,2])) + exp(data[r_num,2,ind]+log(m[2,2])))

        @inbounds data[n_num,1,ind] = a+c
        @inbounds data[n_num,2,ind] = b+d
        res += a+b+c+d

    end # for
    return res
end # function

function expo2(p, t, r)
    return [-p*exp(-t*r) p*exp(-t*r);
            (1 - p)*exp(-t*r) (p - 1)*exp(-t*r)]
end
function expo1(p,t,r)
    return [exp(-t*r)+p*(1-exp(-r*t)) (1-p)*exp(-r*t);
            p*(1-exp(-r*t)) exp(-r*t)+(1-p)*(1-exp(-r*t))]
end

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
