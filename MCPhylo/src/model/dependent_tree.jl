
function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(AbstractString["Tree height", "Tree length"], n_names)
end

function names(d::TreeLogical, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(AbstractString["Tree height", "Tree length"], n_names)
end


function unlist(d::TreeStochastic)
    y = tree_height(d.value)
    x = vec([n.height for n in post_order(d.value) if n.root !== true])
    vcat(y, tree_length(d.value), x)
end

function unlist(d::TreeLogical)
    y = tree_height(d.value)
    x = vec([n.height for n in post_order(d.value) if n.root !== true])
    vcat(y, tree_length(d.value), x)
end


function unlist(s::AbstractStochastic, x::N, transform::Bool=false)::Vector{N}  where N <: AbstractNode
    [s.value]
end
