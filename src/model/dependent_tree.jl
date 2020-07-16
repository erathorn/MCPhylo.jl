
function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(n_names,AbstractString["Tree height", "Tree length"])
end

function names(d::TreeLogical, nodekey::Symbol)
    n_names = [n.num for n in post_order(d.value) if n.root !== true]
    sort!(n_names)
    n_names = vec(AbstractString["node "*string(n) for n in n_names])
    AbstractString["Tree height", "Tree length"]
    vcat(n_names, AbstractString["Tree height", "Tree length"])
end


function unlist(root::Node)
    x = node_height_vec(root)
    vcat(x, tree_length(root))
end

function unlist(d::TreeStochastic)
    unlist(d.value)
end

function unlist(d::TreeLogical)
    unlist(d.value)
end


function unlist(s::AbstractStochastic, x::N, transform::Bool=false)::Vector{N}  where N <: GeneralNode
    [s.value]
end
