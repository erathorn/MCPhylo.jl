module Tree_moves

include("./Tree_Basics.jl")
using Markdown
using ..Tree_Basics: Node, post_order, set_binary!, add_child!, remove_child!

export NNI!

"""
    NNI!(root::Node)

This function performs an inplace nearest neighbour interchange operation on the
tree which is supplied.
"""
function NNI!(root::Node)
    post_order_trav = post_order(root)
    target::Node = Node(1.0, [0.0], Node[], 0, true, 0.0, "0")
    while true

        target = rand(post_order_trav)
        println(target.name)
        if target.nchild != 0
            if target.child[1].nchild !=0
                if target.child[2].nchild !=0
                    break
                end
            end # if
        end # if
    end # end while

    if rand([1,2]) == 1
        child1 = remove_child!(target, 1)
        child2 = remove_child!(target, 2)

        gchild1 = remove_child!(child1, 1)
        gchild2 = remove_child!(child1, 1)
    else
        child1 = remove_child!(target, 2)
        child2 = remove_child!(target, 1)
        gchild1 = remove_child!(child1, 1)
        gchild2 = remove_child!(child1, 1)
    end # if

    add_child!(target, child1)
    add_child!(target, gchild1)
    add_child!(child1, child2)
    add_child!(child1, gchild2)

    set_binary!(root)

end # function


end # module Tree_moves
