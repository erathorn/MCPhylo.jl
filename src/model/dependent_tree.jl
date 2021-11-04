
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
    s.value
end
