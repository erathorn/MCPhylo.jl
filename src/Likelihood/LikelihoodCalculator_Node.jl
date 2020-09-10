
"""
FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64},
                             data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}
    r::Float64 = 1.0
    mu =  1.0 / (2.0 * pi_ * (1-pi_))
    mml = calc_trans.(blv, pi_, mu, r)
    #println(mml[3])
    root_node = last(tree_postorder)
    rns = zeros(Float64, 1, n_c)
    res::Array{Float64, 2} = ones(Float64, 2, n_c)
    #@views
    for node in tree_postorder
        if node.nchild > 0
            node.data = node_loop(node, mml)
            #println(sum(node.data),"  ", node.num)
            #println("--")

            if !node.root
                scaler = Base.maximum(node.data, dims=1)
                rns = rns + scaler
                node.data = node.data .- scaler
            end #if

        end #if
    end # for

    #println(sum(root_node.data))
    #println(sum(rns))
    return MCPhylo.likelihood_root(exp.(root_node.data), pi_, rns)
end # function

@inline function likelihood_root(root::Array{A,2}, pi_::S, rns::Array{A,2})::A where {S<:Real, A<:Real}
    sum(log.(sum(root .* Array([pi_, 1.0-pi_]), dims=1))+rns)
end

function node_loop(node::N, mml::Array{Array{Float64, 2},1})::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    out = zeros(size(node.data))
    #@inbounds @views
    for child in node.children
            out = out .+ log.(exp.(mml[child.num])*exp.(child.data))
    end

    out
end


function node_loop2(node::T, mml::Array{Array{S, 2},1})::Array{S,2} where {T<:Node, S<:Real}
    reduce(pointwise_reduce2, bc2.(node.children, Ref(mml)))
end

function pointwise_reduce2(x::S, y::Array{T})::Array{T} where {S<:Array, T<:Real}
    x.*y
end

function bc2(node::T, mml::Array{Array{S,2},1})::Array{S} where {T<:Node, S<:Real}
    mml[node.num]*node.data
end

@inline function pointwise_reduce(x::Array{T,N}, y::Array{T,N})::Array{T,N} where {N, T<:Real}
    x.*y
end

@inline function bc(node::N, mml::Array{Array{Float64,2},1})::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    mml[node.num]*node.data
end


function calc_trans(time::Float64, pi_::S, mu::Float64, r::Float64)::Array{Float64,2} where {S<:Real}

    v = exp(-r*time*mu)
    v1 = pi_ - (pi_ - 1)*v
    v2 = pi_ - pi_* v
    v3 = (pi_ - 1)*(v - 1)
    v4 = pi_*(v - 1) + 1
    mm::Array{Float64}=[v1 v3;
                        v2 v4]
    return log.(mm)
end
