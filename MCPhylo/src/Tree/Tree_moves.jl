
#TODO: RSPR

"""
    NNI(root::T, target::T, lor::Bool)::Int64   where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The parameter `target` specifies the node which performs the interchange
move using the left or right child of the target node. If the left child should
be used `lor=true`.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T, target::T, lor::Bool)::Int64  where T<:AbstractNode
    # NNI move would be illegal
    if target.nchild === 0 || target.root
        return 0
    end # if

    parent::T = get_mother(target)
    sister::T = get_sister(target)

    ychild::T = remove_child!(target, lor)
    xchild::T = remove_child!(parent, sister)

    add_child!(target, sister)
    add_child!(parent, ychild)

    set_binary!(root)

    return 1

end # function

"""
    NNI!(root::T, target::Int64)::Int64  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The target is identified by the number of the target node.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T, target::Int64)::Int64  where T<:AbstractNode
   tn::T = find_num(root, target)
   lor::Bool = 0.5 > rand()
   NNI!(root, tn, lor)
end #function

"""
    NNI!(root::T)::Int64  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The target is identified by the number of the target node.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T)::Int64  where T<:AbstractNode
    n = rand(1:size(root)[1])
    tn::T = find_num(root, n)
    lor::Bool = 0.5 > rand()
    NNI!(root, tn, lor)
end #function


"""
    slide!(root::Node)

This functin performs a slide move on an intermediate node. The node is moved
upwards or downwards on the path specified by its mother and one of its
daughters.
"""
function slide!(root::Node)
    throw("I need repair")
    target::Node = Node(1.0, [0.0], Node[], 0, true, 0.0, "0")
    while true
        target = random_node(root)
        # check if target is not a leave and that its grand daughters are also
        # no leaves
        if target.nchild != 0
            if target.child[1].nchild !=0
                if target.child[2].nchild !=0
                    break
                end
            end # if
        end # if
    end # end while

    # proportion of slide move is randomly selected
    proportion::Float64 = rand(Uniform(0,1))

    # pick a random child
    child::Node = target.child[rand([1,2])]

    # calculate and set new values
    move!(target, child, proportion)

end # function slide!

"""
    swing!(root::Node)

This function performs a swing node. A random non-leave node is selected and
moved along the path specified by its two children.
"""
function swing!(root::Node)
    throw("I need repair")
    target::Node = Node(1.0, [0.0], Node[], 0, true, 0.0, "0")
    while true
        target = random_node(root)
        # check if target is not a leave
        if target.nchild != 0
            break
        end # if
    end # end while

    proportion::Float64 = rand(Uniform(0,1))

    child1 = target.child[1]
    child2 = target.child[2]

    # calculate and set new values
    move!(child1, child2, proportion)
end # function swing!


"""
    randomize!(root::Node, num::Int64=100)::nothing

This function randomizes the tree topology by performing a number of nearest
neighbour interchange (NNI) moves. The number of NNI moves is specified in
the parameter num.
"""
function randomize!(root::T, num::Int64=100)::Nothing where T <:AbstractNode
    n_nodes = size(root)[1]
    i = 0
    while i < num
        n = rand(1:n_nodes)
        NNI!(root, n)
        i += 1
    end

end



"""
    move!(node1::Node, node2::Node, proportion::Float64)

Change the incomming length of node1 and node2 while keeping there combined length
constant.
"""
function slide!(node1::T, node2::T, proportion::Float64) where T <:AbstractNode
    total::Float64 = node1.inc_length + node2.inc_length
    fp::Float64 = total*proportion
    sp::Float64 = total-fp
    node1.inc_length = fp
    node2.inc_length = sp
end # function slide!
