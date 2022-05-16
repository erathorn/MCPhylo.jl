
function update!(d::Stochastic{<:GeneralNode}, m::Model)
    d.distr = d.eval(m)
    d
end

function update!(d::Logical{<:GeneralNode}, m::Model)
    d.value = d.eval(m)
    d
end

function names(d::Stochastic{<:GeneralNode}, nodekey::Symbol)
    AbstractString["Tree_length["*string(nodekey)*"]"]
end

function names(d::Logical{<:GeneralNode}, nodekey::Symbol)
    AbstractString["Tree_length["*string(nodekey)*"]"]
end

function unlist(root::N) where N<:GeneralNode
    [tree_length(root)]
end

function unlist(d::T) where T <: TreeVariate
    unlist(d.value)
end

function unlist_tree(root::N) where N<:GeneralNode
    [tree_length(root)]
end


function unlist_tree(s::TreeVariate)
    unlist_tree(s.value)
end

function unlist(s::TreeVariate, x::N, transform::Bool=false)  where N <: GeneralNode
    transform ? link_sub(s.distr, s.value) : s.value
end


function link_sub(d::TreeDistribution, X::N) where N<:GeneralNode
    com = to_covariance(X)
    another_unlist(Bijectors.pd_link(com))
end

function relistlength_sub(d::TreeDistribution, s::Stochastic{<:AbstractArray{Float64, N} where N}, X::AbstractArray{<:Real})
    another_relist(X)
end

function relistlength(d::TreeDistribution, s::Stochastic{<:AbstractArray{Float64, N} where N}, X::AbstractArray{<:Real} )
    another_relist(X)
end

function relistlength(d::TreeDistribution, X::AbstractArray{<:GeneralNode})
    relistlength(d, X[1])
end

function relistlength(d::TreeDistribution, X::AbstractArray{<:Real})
    another_relist(X)
end

function invlink_sub(d::TreeDistribution, X::AbstractArray)
    Y = Bijectors.replace_diag(exp, X)
    Y = Bijectors.getpd(Y)
    #Y[Y .< 0] .= 0
    t = cov2tree(Y, string.(1:size(Y,1)),collect(1:size(Y,1)))
    set_binary!(t)
    number_nodes!(t)
    blv = get_branchlength_vector(t)
    set_branchlength_vector!(t, molifier.(blv, 0.003))
    t
end

function another_unlist(X::AbstractArray)
    n = size(X, 1)
    y = [X[i, j] for i = 1:n for j in 1:n if i >= j]
    y
end

function another_relist(X::AbstractArray)
    d = length(X)
    n = 0.5*(sqrt(8*d+1)-1)
    Y = [i >= j ? X[getix(i, j)] : X[getix(j, i)] for j = 1:n, i = 1:n]
    k = Int(n * (n + 1) / 2)
    (Y, k)
end



