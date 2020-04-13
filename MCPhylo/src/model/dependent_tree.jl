
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


function unlist(root::Node)
    y = tree_height(root)
    x = node_height_vec(root)[1:end-1]
    # tester = [n.height for n in get_leaves(root)]
    # if any(tester .!= 0.0)
    #     println("never happen")
    #     println(tester)
    # end
    # x = vec([n.height for n in post_order(root) if n.root != true])
    vcat(y, tree_length(root), x)
end

function unlist(d::TreeStochastic)
    unlist(d.value)
end

function unlist(d::TreeLogical)
    unlist(d.value)
end


function unlist(s::AbstractStochastic, x::N, transform::Bool=false)::Vector{N}  where N <: AbstractNode
    [s.value]
end
