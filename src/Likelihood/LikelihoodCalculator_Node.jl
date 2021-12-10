
"""
    FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
The function is written such that it is differentiable by Zygote 0.5.3.

Returns log-likelihood as a Real number.

* `tree_postorder` : Vector of Nodes.

* `pi` : Real number used in calculation.

* `rates` : Vector of Floats, used in calculation.

* `data` : Array of Floats, currently a placeholder.

* `n_c` : Int64, nchar value derived from tree and dataset.

* `blv` : Vector of Floats, branchlength vector derived from tree.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Float64,
                             data::Array{Float64,3}, n_c::Int64, blv::Array{Float64,1})::Float64 where {N<:GeneralNode}

    mu::Float64 =  1.0 / (2.0 * prod(pi_))
    mutationArray::Vector{Array{Float64, 2}} = calc_trans.(blv, pi_[1], mu, rates)
    root_node::N = last(tree_postorder)
    rns::Array{Float64} = zeros(Float64, 1, n_c)

    @views for node in tree_postorder
        if node.nchild > 0
            res = node_loop(node, mutationArray, data)
            if !node.root
                scaler::Array{Float64} = maximum(res, dims=1)
                rns = rns .+ log.(scaler)
                res = res ./ scaler
            end #if
            
            data[:,:,node.num] = res
        end #if
    end # for
    return likelihood_root(data[:,:,root_node.num], pi_, rns)
end # function

@inline function likelihood_root(root::Array{Float64,2}, pi_::Array{Float64}, rns::Array{Float64,2})::Float64
    sum(log.(sum(root .* pi_, dims=1))+rns)
end

function node_loop(node::N, mutationArray::Vector{Array{Float64, 2}}, data)::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    # creating a new array is necessary because Zygote can not differentiate otherwise
    out = ones(size(data[:,:,node.num]))
    @inbounds @views for child in node.children
            out = out .* BLAS.gemm('N','N',mutationArray[child.num], data[:,:,child.num])
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
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Float64,
                      U::A, D::V, Uinv::A, mu::Float64,
                     data::Array{Float64,3}, substitutionModel::Function,  c_grad::Bool = true) where {N<:GeneralNode, M<:Number, A<:AbstractArray, V<:AbstractVector}
    Nbases, Nsites, Nnodes = size(data)
    
    
    grv::Vector{Float64} = Vector{Float64}(undef, Nnodes-1)
    Down::Array{Float64,3} = similar(data)
    pre_order_partial::Array{Float64,3} = similar(data)
    ll = fels_ll(tree_postorder, data, D, U, Uinv, rates, mu, Down, pi_, substitutionModel)
    
    if c_grad
        grv = fels_grad(tree_postorder, data, D, U, Uinv, rates, mu,
                        Down, pi_, pre_order_partial, grv, substitutionModel)
    end

    return ll, grv, data, Down, pre_order_partial

end

function fels_grad(tree_postorder::Vector{N}, data::Array{Float64,3},
         D::V, U::A, Uinv::A, rate::Float64, mu::Float64, Down::Array{Float64,3},
         pi_::Array{Float64}, pre_order_partial::Array{Float64,3}, grv::Vector{Float64}, substitutionModel::Function)::Vector{Float64} where {N <: GeneralNode, M <: Number, 
                                                                                                                            A<:AbstractArray, V<:AbstractVector}

    root_node = last(tree_postorder)
    pre_order_partial[:, :, root_node.num] .= pi_

    @inbounds @views for node in reverse(tree_postorder)[2:end]
           mother = get_mother(node)
            pre_order_partial[:, :, node.num] .= pre_order_partial[:, :, mother.num]
            @simd for c in 1:mother.nchild
                child = mother.children[c]
                if child.num != node.num
                    pre_order_partial[:, :, node.num] .*= Down[:, :, child.num]
                end
            end
                
            ptg = U * ((D .* (rate * mu)) .* diagm(exp.(D .* (rate*mu*node.inc_length)))) * Uinv
            gradi = sum(pre_order_partial[:, :,  node.num] .* (ptg * data[:, :, node.num]),dims=1)
            tr = calculate_transition(substitutionModel, rate, mu, node.inc_length, U, Uinv, D, pi_)
            pre_order_partial[:, :,  node.num] .= transpose(tr) * pre_order_partial[:, :, node.num]
            grv[node.num] = sum(gradi ./ sum(pre_order_partial[:, :,  node.num] .* data[:, :,  node.num], dims=1))
            
            if node.nchild > 0
                pre_order_partial[:, :,  node.num] ./= maximum(pre_order_partial[:, :,  node.num], dims = 1)
            end
    end # for
    grv
end

function fels_ll(tree_postorder::Vector{N}, data::Array{Float64,3},
          D::V, U::A, Uinv::A,
         rate::Float64, mu::Float64,Down::Array{Float64,3},
         pi_::Array{Float64}, substitutionModel::Function)::Float64 where {N <: GeneralNode, M<:Number, A<:AbstractArray, V<:AbstractVector}

    ll::Float64 = 0.0
    @inbounds @views for node in tree_postorder
        if node.nchild > 0
            data[:, :, node.num] .= 1.0
            for child in node.children
                mul!(Down[:, :, child.num], calculate_transition(substitutionModel, rate, mu, child.inc_length, U, Uinv, D, pi_), data[:, :, child.num])
                data[:, :, node.num] .*= Down[:, :, child.num]
            end
            
            if !node.root   
                ll += sum(log.(maximum(data[:, :, node.num], dims=1)))
                data[:, :, node.num] ./= maximum(data[:, :, node.num], dims=1)
                #ll += sum(log.(scaler))
            else
                ll += sum(log.(sum(data[:, :, last(tree_postorder).num] .* pi_, dims=1)))
            end
        end #if
    end # for
    ll
end


function recursive_ll(node::N, data::Array{Float64,3},
                        D::V, U::A, Uinv::A,
                    rate::Float64, mu::Float64,Down::Array{Float64,3},
                    pi_::Array{Float64}, substitutionModel::Function)::Float64 where {N <: GeneralNode, M<:Number, A<:AbstractArray, V<:AbstractVector}
    ll = 0.0
    if node.nchild > 0
        data[:, :, node.num] .= 1.0
        for child in node.children
            ll += recursive_ll(child, data, D, U, Uinv, rate, mu, Down, pi_, substitutionModel)
   
            tr = calculate_transition(substitutionModel, rate, mu, child.inc_length, U, Uinv, D, pi_)
            Down[:, :, child.num] .= tr * data[:, :, child.num]            
            data[:, :, node.num] .*= Down[:, :, child.num]
        end
        
        if !node.root   
            scaler = maximum(data[:, :, node.num], dims=1)
            data[:, :, node.num] ./= scaler
            ll += sum(log.(scaler))
        else
            ll += sum(log.(sum(data[:, :, node.num] .* pi_, dims=1)))
        end
    end 
    return ll
end