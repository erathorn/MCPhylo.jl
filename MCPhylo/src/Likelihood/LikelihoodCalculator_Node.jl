
@views function Felsenstein_Recursion(root::N, pi_::T, rates::Vector{Float64}, data::Array{S,3},
                               n_c::Int64, trmats::Array{Array{S, 2},1},
                               od_scaler::Array{S, 2})::Tuple{Array{S,2},Array{S,2}}  where {S<:Real, T<:Real, N<:AbstractNode}


   if root.nchild == 0
       return data[:,:,root.num], od_scaler
   end

   # here happens the intresting stuff
   child_data::Array{Tuple{Array{S,2},Array{S,2}},1} = Felsenstein_Recursion.(root.children, Ref(pi_), Ref(rates), Ref(data), n_c, Ref(trmats), Ref(od_scaler))

   num_v::Array{Array{S,2},1} = [trmats[node.num] for node in root.children]
   rv::Array{S,2} = reduce(pointwise_reduce, num_v .* first.(child_data))
   od_scaler = sum(last.(child_data))
   if !root.root
       scaler = Base.maximum(rv, dims=1)
       od_scaler = od_scaler + log.(scaler)
       rv = rv ./ scaler
   end
   return rv, od_scaler

end


function FelsensteinFunction(root::N, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Array{S})::S where {S<:Real, T<:Real, N<:AbstractNode}
    r::Float64 = 1.0
    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = calc_trans.(blv, pi_, mu, r)
    od_scaler = zeros(S, 1, n_c)

    rv, od_scaler = Felsenstein_Recursion(root, pi_, rates, data, n_c, mml, od_scaler)


    sum(log.(sum(rv .* Array([pi_, 1.0-pi_]), dims=1)) + od_scaler)
end

"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64)::Float64

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64},
                             data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:Node{<:Real,Float64,Float64,<:Integer}}
    r::Float64 = 1.0
    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = calc_trans.(blv, pi_, mu, r)

    root_node = last(tree_postorder)
    root_node.scaler = zeros(1, n_c)
    #rns = zeros(S, 1, n_c)
    @views for node in tree_postorder
        if node.nchild > 0
            node.data = node_loop(node, mml)
            if !node.root
                node.scaler = Base.maximum(node.data, dims=1)
                root_node.scaler = root_node.scaler + log.(node.scaler)
                node.data = node.data ./ node.scaler
            end #if
        end #if
    end # for

    return likelihood_root(root_node.data, pi_, root_node.scaler)#res
end # function

@inline function likelihood_root(root::Array{A,2}, pi_::S, rns::Array{A,2})::A where {S<:Real, A<:Real}
    sum(log.(sum(root .* Array([pi_, 1.0-pi_]), dims=1))+rns)
end

function node_loop(node::N, mml::Array{Array{Float64, 2},1})::Array{Float64,2} where {N<:Node{<:Real,Float64,Float64,<:Integer}, S<:Real}
    r = bc.(node.children, Ref(mml))
    reduce(pointwise_reduce, r)
end

@inline function pointwise_reduce(x::Array{T,N}, y::Array{T,N})::Array{T,N} where {N, T<:Real}
    x.*y
end

@inline function bc(node::N, mml::Array{Array{Float64,2},1})::Array{Float64,2} where {N<:Node{<:Real,Float64,Float64,<:Integer}, S<:Real}
    r = node.data
    s = mml[node.num]
    s*r
end


function calc_trans(time::Float64, pi_::S, mu::Float64, r::Float64)::Array{Float64,2} where {S<:Real}

    v = exp(-r*time*mu)
    v1 = pi_ - (pi_ - 1)*v
    v2 = pi_ - pi_* v
    v3 = (pi_ - 1)*(v - 1)
    v4 = pi_*(v - 1) + 1
    mm::Array{Float64}=[v1 v3;
                  v2 v4]
    return mm
end

function nf(node::T, blv::Vector, pi_::S, mu::Float64, r::Float64, mm::Array{Float64,2}) where {T<:AbstractNode, S<:Number}
    md::Array{Float64, 2} = node.data
    @inbounds m_num::Int64 = node.num
    @inbounds minc::Float64 = blv[m_num]
    v::Float64 = exp(-r*minc*mu)
    @inbounds mm[1] = pi_ - (pi_ - 1)*v
    @inbounds mm[2] = pi_ - pi_* v
    @inbounds mm[3] = (pi_ - 1)*(v - 1)
    @inbounds mm[4] = pi_*(v - 1) + 1

    mm*md
end

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
"""
function CondLikeInternal(pi_::T, l_data::Array{Float64}, r_data::Array{Float64}, linc::Float64, rinc::Float64)::Array where T<:Number#::Float64


        r::Float64 = 1.0
        mu::Float64 =  1.0 / (2.0 * pi_ * (1-pi_))
        v::Float64 = exp(-r*linc*mu)
        m11::Float64 = pi_ - (pi_ - 1)*v
        m22::Float64 = pi_*(v - 1) + 1
        m21::Float64 = pi_ - pi_* v
        m12::Float64 = (pi_ - 1)*(v - 1)
        mm = [m11 m12; m21 m22]

        v1::Float64 = exp(-r*rinc*mu)
        m11b::Float64 = pi_ - (pi_ - 1)*v1
        m22b::Float64 = pi_*(v1 - 1) + 1
        m21b::Float64 = pi_ - pi_* v1
        m12b::Float64 = (pi_ - 1)*(v1 - 1)
        mm2 = [m11b m12b; m21b m22b]

        return (mm*l_data) .* (mm2*r_data)


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
