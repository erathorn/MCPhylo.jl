
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
    if target.nchild == 0 || target.root
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


### Experimental rerooting



function reroot(root::T, new_root::String)::T where T<:Node

    new_tree = deepcopy(root)
    root_node = find_by_name(new_tree, new_root)

    mother = root_node.mother

    recursive_invert(mother, root_node)

    root_node.root = true
    new_tree.root = false


    set_binary!(root_node)
    number_nodes!(root_node)
    return root_node
end


function recursive_invert(old_mother::T, old_daughter::T)::T where T
    if old_mother.root == true
        # arrived at the root
        od = remove_child!(old_mother, old_daughter)
        add_child!(od, old_mother)
        return od
    end
        od1 = recursive_invert(old_mother.mother, old_mother)
        od = remove_child!(od1, old_daughter)
        add_child!(od, od1)
        return od
end

"""
    SPR(original_root::Node)::AbstractNode
Performs SPR on tree; takes a copy of root of the tree;
Returns a copy of root of altered tree
"""


function SPR(original_root::Node)
    root = deepcopy(original_root)
    SPR!(root)
    return root
end #func

"""
        SPR!(root::Node)::AbstractNode
    Performs SPR on tree in place; takes reference to root of tree, boolean value necessary to determine if tree should be treated as binary or not
    Returns reference to root of altered tree
"""

function SPR!(root::Node)
    if length(post_order(root)) <= 2
        error("The tree is too small for SPR")
    end #if
    binary = check_binary(root)
    spr_tree = binary ? perform_spr(root) : throw("Not yet implemented for not binary trees")
    return spr_tree
end #function

"""
        risky_SPR(root::Node)::AbstractNode
    Performs SPR on tree in place; takes reference to root of tree, boolean value necessary to determine if tree should be treated as binary or not
    Returns copy of root of altered tree. Does not check for correct formatting of tree.
"""
function risky_SPR(original_root::Node)
    root = deepcopy(original_root)
    return risky_SPR!(root)
end #function


"""
        risky_SPR!(root::Node)::AbstractNode
    Performs SPR on tree in place; takes reference to root of tree, boolean value necessary to determine if tree should be treated as binary or not
    Returns reference to root of altered tree. Does not check for correct formatting of tree.
"""
function risky_SPR!(root::Node)
    return perform_spr(root)
end #func

"""
    perform_spr_binary(root::Node)::AbstractNode
performs SPR on binary tree
"""
function perform_spr(root::Node)
    # find node to move
    available = [n.num for n in post_order(root)]
    n = rand(available)
    tn::Node = MCPhylo.find_num(root, n) #his is the root of the subtree which will be moved
    while tn.root || tn.mother.root
        n = rand(available)
        tn = MCPhylo.find_num(root, n) #his is the root of the subtree which will be moved
    end # while
    tn_mother = tn.mother
    tn_sister = MCPhylo.get_sister(tn)
    tn_gm = tn_mother.mother
    remove_child!(tn_gm, tn_mother)
    remove_child!(tn_mother, tn_sister)
    add_child!(tn_gm, tn_sister)
    # find target
    available = [n.num for n in post_order(root)]
    n = rand(available)
    target::Node = MCPhylo.find_num(root, n) #this is the target of the movement
    while target.root
        n = rand(available)
        target = MCPhylo.find_num(root, n)
    end # while
    target_mother = target.mother
    remove_child!(target_mother, target)
    add_child!(target_mother, tn_mother)
    add_child!(tn_mother, target)
    return root
end #func
