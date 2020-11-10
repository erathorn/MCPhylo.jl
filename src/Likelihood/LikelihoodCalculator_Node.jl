
"""
FelsensteinFunction(tree_postorder::Vector{N}, pi_::T, rates::Vector{Float64}, data::Array{Float64,3}, n_c::Int64, blv::Vector{Float64}) where {T<:Real, N<:GeneralNode}

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.

The function is written such that it is differentiable by Zygote 0.5.3.
"""
function FelsensteinFunction(tree_postorder::Vector{N}, pi_::Array{Float64}, rates::Float64,
                             data::Array{Float64,3}, n_c::Int64, blv::Array{Float64,1})::Float64 where {N<:GeneralNode}

    mu::Float64 =  1.0 / (2.0 * prod(pi_))
    mutationArray::Array{Float64} = reshape(vcat(calc_trans.(blv, Ref(pi_), mu, rates)...) ,
                            (length(blv), length(pi_), length(pi_)))
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

function node_loop(node::N, mutationArray::Array{Float64, 3})::Array{Float64,2} where {N<:GeneralNode, S<:Real}
    # creating a new array is necessary because Zygote can not differentiate otherwise
    #out = ones(size(node.data))
    #@inbounds @views for child in node.children
    #        out = out .* (mutationArray[child.num]*child.data)
    #end
    #out
    mapreduce(y->mutationArray[y.num, :, :]*y.data, (x,z)->x .* z, node.children)
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

function calc_trans(time::Float64, pi_::S, mu::Float64, r::Array{T,1})::Array{Float64,3} where {S<:Real, T<:Real}
    permutedims(reshape(hcat(map(i-> MCPhylo.calc_trans(time, pi_, mu, i), r)...), 2,2,length(r)), [3,1,2])
end
