
"""
FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.

The function is written such that it is differentiable by Zygote 0.5.3.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Float64,
                             data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = calc_trans.(blv, pi_, mu, rates)

    root_node = last(tree_postorder)
    rns = zeros(Float64, 1, n_c)

    @views for node in tree_postorder
        if node.nchild > 0
            node.data = node_loop(node, mml)
            if !node.root
                scaler = maximum(node.data, dims=1)
                rns = rns + log.(scaler)
                node.data = node.data ./ scaler
            end #if
        end #if
    end # for

    return likelihood_root(root_node.data, pi_, rns)
end # function

@inline function likelihood_root(root::Array{A,2}, pi_::S, rns::Array{A,2})::A where {S<:Real, A<:Real}
    sum(log.(sum(root .* Array([pi_, 1.0-pi_]), dims=1))+rns)
end

function node_loop(node::N, mml::Array{Array{Float64, 2},1})::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    # creating a new array is necessary because Zygote can not differentiate otherwise
    out = ones(size(node.data))
    @inbounds @views for child in node.children
            out = out .* (mml[child.num]*child.data)
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
