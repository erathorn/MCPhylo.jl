
function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function update!(d::TreeLogical, m::Model)
    d.value = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    AbstractString["Tree_height["*string(nodekey)*"]" , "Tree_length["*string(nodekey)*"]"]
end

function names(d::TreeLogical, nodekey::Symbol)
    AbstractString["Tree_height["*string(nodekey)*"]" , "Tree_length["*string(nodekey)*"]"]
end


function unlist(root::N) where N<:GeneralNode
    [tree_height(root), tree_length(root)]
end


function unlist(d::T) where T <: TreeVariate
    unlist(d.value)
end


function unlist_tree(root::N) where N<:GeneralNode
    [tree_height(root), tree_length(root)]
end


function unlist_tree(s::AbstractTreeStochastic)
    unlist_tree(s.value)
end

function unlist(s::AbstractTreeStochastic, x::N, transform::Bool=false)  where N <: GeneralNode
    s.value
end
