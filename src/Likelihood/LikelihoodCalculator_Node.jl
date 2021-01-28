
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

"""
    FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Array{Float64},
                         U::Array{Float64,2}, D::Array{Float64}, Uinv::Array{Float64,2},
                         data::Array{Float64,4}, c_grad::Bool = true) where {N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm. If `c_grad` equals `true` (default) the analytic gradient
regarding the branch lengths of the tree gets computed as well.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Array{Float64},
                     U::Array{Float64,2}, D::Array{Float64}, Uinv::Array{Float64,2}, mu::Float64,
                     data::Array{Float64,4}, c_grad::Bool = true) where {N<:GeneralNode}
    Nbases, Nsites, Nrates, Nnodes = size(data)
    mutationArray::Array{Float64,4} = Array{Float64,4}(undef, Nbases, Nbases, Nrates, Nnodes-1)
    grv::Vector{Float64} = Vector{Float64}(undef, Nnodes-1)
    Down::Array{Float64,4} = similar(data)
    ll = fels_ll(tree_postorder, data, D, U, Uinv, rates, mu, Nrates, Nsites, Down, pi_, mutationArray)
    if c_grad
        grv = fels_grad(tree_postorder, data, D, U, Uinv, rates, mu, Nrates, Nsites, Nnodes, Down, pi_, mutationArray)
    end

    return ll, grv

end

function fels_grad(tree_postorder::Vector{N}, data::Array{Float64,4},
         D::Array{Float64,1}, U::Array{Float64,2}, Uinv::Array{Float64,2},
         rates::Array{Float64,1}, mu::Float64, Nrates::Int64, Nsites::Int64, Nnodes::Int64, Down::Array{Float64,4},
         pi_::Array{Float64}, mutationArray::Array{Float64,4})::Vector{Float64} where {N <: GeneralNode}

    root_node::N = last(tree_postorder)
    pre_order_partial::Array{Float64,4} = similar(data)
    pre_order_partial[:, :, :, root_node.num] .= pi_
    scaler::Array{Float64, 2} = Array{Float64,2}(undef, 1, Nsites)
    gradi::Array{Float64, 2} = Array{Float64,2}(undef, 1, Nsites)
    grv::Vector{Float64} = Vector{Float64}(undef, Nnodes-1)
    ptg::Array{Float64,2} = similar(mutationArray[:, :, 1, 1])

    @inbounds @views for node in reverse(tree_postorder)[2:end]
           mother::N = node.mother
           @inbounds @views for r in 1:Nrates
               pre_order_partial[:, :, r, node.num] .= pre_order_partial[:, :, r, mother.num]
               @inbounds @views for child in mother.children
                   if child.num != node.num
                       pre_order_partial[:, :, r, node.num] .*= Down[:, :, r, child.num]
                   end
               end
           end
           tg::Float64 = 0.0
           @inbounds @views for r in 1:Nrates
               BLAS.gemm!('N', 'N', 1.0, BLAS.symm('R', 'L', (D .* (rates[r] * mu)) .* diagm(exp.(D .* (rates[r]*mu*node.inc_length))), U), Uinv, 0.0, ptg)
               gradi = sum(pre_order_partial[:, :, r, node.num] .* BLAS.gemm('N', 'N', ptg , data[:, :, r, node.num]), dims=1)
               pre_order_partial[:, :, r, node.num] .= BLAS.gemm('T', 'N', mutationArray[:, :, r, node.num] , pre_order_partial[:, :, r, node.num])
               tg += sum(gradi ./ sum(pre_order_partial[:, :, r, node.num] .* data[:, :, r, node.num], dims=1))
               if node.nchild > 0
                   scaler = maximum(pre_order_partial[:, :, r, node.num], dims = 1)
                   pre_order_partial[:, :, r, node.num] ./= scaler
               end
           end
           grv[node.num] = tg
    end # for
    grv
end

function fels_ll(tree_postorder::Vector{N}, data::Array{Float64,4},
         D::Array{Float64,1}, U::Array{Float64,2}, Uinv::Array{Float64,2},
         rates::Array{Float64,1}, mu::Float64, Nrates::Int64, Nsites::Int64, Down::Array{Float64,4},
         pi_::Array{Float64}, mutationArray::Array{Float64,4})::Float64 where {N <: GeneralNode}

    scaler::Array{Float64, 2} = Array{Float64,2}(undef, 1, Nsites)
    ll::Float64 = 0.0
    root_node::N = last(tree_postorder)
    @inbounds @views for node in tree_postorder
        if node.nchild > 0
            data[:, :, :, node.num] .= 1.0
            for child in node.children
                @inbounds @views for r in 1:Nrates
                    BLAS.gemm!('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* (rates[r]*mu*child.inc_length))), U), Uinv, 0.0, mutationArray[:, :, r, child.num])
                    BLAS.gemm!('N','N', 1.0, mutationArray[:, :, r, child.num], data[:, :, r, child.num], 0.0, Down[:, :, r, child.num])
                    data[:, :, r, node.num] .*= Down[:, :, r, child.num]

                end

            end

            if !node.root
                @inbounds @views for r in 1:Nrates
                    scaler .= maximum(data[:, :, r, node.num], dims=1)
                    data[:, :, r, node.num] ./= scaler
                    ll += sum(log.(scaler))
                end
            end
        end #if
    end # for

    @inbounds @views @simd for r in 1:Nrates
        ll += sum(log.(sum(data[:, :, r, root_node.num] .* pi_, dims=1)))
    end
    ll
end
