
"""
FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.

The function is written such that it is differentiable by Zygote 0.5.3.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Float64,
                             data::Array{Float64,3}, n_c::Int64, blv::Array{Float64,1})::Float64 where {N<:GeneralNode}

    mu::Float64 =  1.0 / (2.0 * prod(pi_))
    mutationArray::Vector{Array{Float64, 2}} = calc_trans.(blv, pi_[1], mu, rates)
    root_node::N = last(tree_postorder)
    rns::Array{Float64} = zeros(Float64, 1, n_c)

    @views for node in tree_postorder
        if node.nchild > 0
            res = node_loop(node, mutationArray)
            if !node.root
                scaler::Array{Float64} = maximum(res, dims=1)
                rns = rns .+ log.(scaler)
                res = res ./ scaler
            end #if
            node.data = res
        end #if
    end # for

    return likelihood_root(root_node.data, pi_, rns)
end # function

@inline function likelihood_root(root::Array{Float64,2}, pi_::Array{Float64}, rns::Array{Float64,2})::Float64
    sum(log.(sum(root .* pi_, dims=1))+rns)
end

function node_loop(node::N, mutationArray::Vector{Array{Float64, 2}})::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    # creating a new array is necessary because Zygote can not differentiate otherwise
    out = ones(size(node.data))
    @inbounds @views for child in node.children
            out = out .* BLAS.gemm('N','N',mutationArray[child.num], child.data)
    end
    out
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


function calc_trans(time::Float64, pi_::Array{Float64}, mu::Float64, r::Float64)::Array{Float64,2}
    res::Array{Float64} = repeat(reshape(pi_, (1,2)), 2)
    di = Diagonal(ones(2))
    v = exp(-r*time*mu)
    res = res .*(1.0-v)
    res .+ (di .* v)
end

function Felsenstein_grad(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Float64,
                             data::Array{Float64,3}, n_c::Int64) where {N<:GeneralNode}

     mu::Float64 =  1.0 / (2.0 * prod(pi_))
     mutationArray::Array{Float64,3} = Array{Float64,3}(undef, 2,2,length(tree_postorder)-1)
     root_node::N = last(tree_postorder)
     ll::Float64 = 0.0
     Down = similar(data)

     @inbounds @views for node in tree_postorder
         if node.nchild > 0
             out = ones(size(node.data))
             @inbounds @views for child in node.children
                 mutationArray[:, :, child.num] = calc_trans(child.inc_length, pi_[1], mu, rates)
                 Down[:, :, child.num] = BLAS.gemm('N','N',mutationArray[:, :, child.num], data[:, :, child.num])
                 out .*= Down[:, :, child.num]
             end
             if !node.root
                 scaler::Array{Float64} = maximum(out, dims=1)
                 out ./= scaler
                 ll += sum(log.(scaler))
             end
             #node.data = out
             data[:, :, node.num] = out
         end #if
     end # for
     ll += sum(log.(sum(data[:, :, root_node.num] .* pi_, dims=1)))

     pre_order_partial = similar(data)
     pre_order_partial[:, :, root_node.num] .= pi_

     grv::Array{Float64,1} = Array{Float64,1}(undef, length(tree_postorder)-1)

     Q = ones(2,2)
     Q[diagind(2,2)] .= -1
     Q .*= pi_

     @inbounds @views for node in reverse(tree_postorder)
             if !node.root
                 mother = node.mother
                 pre_order_partial[:, :, node.num] .= pre_order_partial[:, :, mother.num]
                 for child in mother.children
                     if child.num != node.num
                         pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
                     end
                 end

                 ptg::Array{Float64,2} = (rates*mu*Q)*mutationArray[:, :, node.num]
                 gradi::Array{Float64,2} = sum(pre_order_partial[:, :, node.num] .* BLAS.gemm('N', 'N', ptg , data[:, :, node.num]), dims=1)
                 #pre_order_partial[:, :, node.num] .= BLAS.gemm('T', 'N', mutationArray[:, :, node.num] , pre_order_partial[:, :, node.num])
                 BLAS.gemm!('T', 'N',1.0, mutationArray[:, :, node.num] , pre_order_partial[:, :, node.num], 0.0, pre_order_partial[:, :, node.num])
                 grv[node.num] = sum(gradi ./ sum(pre_order_partial[:, :, node.num] .* data[:, :, node.num], dims=1))
                 if node.nchild > 0
                     scaler = maximum(pre_order_partial[:, :, node.num], dims = 1)
                     pre_order_partial[:, :, node.num] ./= scaler
                 end
             end #if

     end # for
     return ll, grv

end
