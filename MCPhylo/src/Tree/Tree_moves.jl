
#TODO: RSPR
"""
    NNI!(root::Node)

This function performs an inplace nearest neighbour interchange operation on the
tree which is supplied.
"""
function NNI!(root::T) where T<:Node

    target::Node = Node("node_name", zeros(Float64, (2, 1)),missing, missing, missing, 0, true, 0.0, "0", 0)
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

end # function NNI!


"""
    NNI(root::Node, target::Node)::nothing

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The parameter `target` specifies the node which performs the interchange
move with a randomly selected neighbour.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T, target::T)::Int64 where T <: Node
    # NNI move would be illegal
    if target.nchild == 0 || target.root
        return 0
    end # if


    parent::Node = get_mother(target)
    sister = get_sister(target)

    rand([1,2]) == 1 ? ychild = remove_child!(target, true) : ychild = remove_child!(target, false)

    xchild = remove_child!(parent, sister)


    add_child!(target, xchild)
    add_child!(parent, ychild)


    set_binary!(root)

    return 1

end # function



"""
    NNI(root::Node, target::Node)::nothing

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The parameter `target` specifies the node which performs the interchange
move with a randomly selected neighbour.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T, target::T, lor::Bool)::Int64  where T<:Node
    # NNI move would be illegal
    if target.nchild === 0 || target.root
        return 0
    end # if


    parent = get_mother(target)

    sister = get_sister(target)

    ychild = remove_child!(target, lor)

    xchild = remove_child!(parent, sister)
    

    add_child!(target, sister)
    add_child!(parent, ychild)


    set_binary!(root)

    return 1

end # function


function NNI!(root::T, target::Int64, lor::Bool)::Int64  where T<:Node
   tn = find_num(root, target)
   NNI!(root, tn, lor)
end #function


function randomize!(root::Node, num::Int64=100)
    nnodes = size(root)[1]
    i = 0
    while i < num
        n = rand(1:nnodes)
        lor = 0.5 > rand()
        NNI!(root, n, lor)
        i+=1
    end
end

"""
    slide!(root::Node)

This functin performs a slide move on an intermediate node. The node is moved
upwards or downwards on the path specified by its mother and one of its
daughters.
"""
function slide!(root::Node)
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
    NNI(mat::Array{Float64,2})::nothing

documentation
"""
function NNI!(mat::Array{Float64,2})::Array{Float64,2}
    leaves::Array{Float64,1} = get_leaves(mat)
    l::Int64 = size(mat)[1]

    target::Int64 = 0

    while true
        target = rand(collect(1:l))
        if !(target in leaves)
            if length(intersect!(vec(get_neighbours(mat[target,:])), vec(leaves))) == 0
                break
            end
        end
    end

    ac = []
    ch::Array{Int64} = get_neighbours(mat[target,:])
    for c in 1:2
        push!(ac, get_neighbours(mat[ch[c],:]))
    end

    if 0.5 > rand()
        # swap ac[1,1] with ac[2, 1]
        swap_cols(mat, ac[1][1], ac[2][1])
    else
        # swap ac[1,2] with ac[2, 1]
        swap_cols(mat, ac[1][2], ac[2][1])
    end

    return mat

end # function


function swap_cols(mat::Array{Float64, 2}, ind::Int64, jnd::Int64)
    l::Int64 = size(mat)[1]
    @assert (1 <= ind <= l) && (1 <= jnd <= l)
    @inbounds for i in 1:l
        mat[i,ind], mat[i,jnd] = mat[i,jnd], mat[i,ind]
    end
    return mat
end


"""
    NNI(mat::Array{Float64,2})::nothing

documentation
"""
function NNI!(mat::Array{Float64,2}, target::Int64)::Int64
    leaves::Vector{Int64} = get_leaves(mat)
    l::Int64 = size(mat)[1]

    if (target in leaves)|(length(intersect!(get_neighbours(mat[target,:]), vec(leaves))) != 0)

        return 0
        #end

    end


    ac = []
    ch::Array{Int64} = get_neighbours(mat[target,:])

    for c in 1:2
        push!(ac, get_neighbours(mat[ch[c],:]))
    end

    if rand([true,false])
        # swap ac[1,1] with ac[2, 1]
        swap_cols(mat, ac[1][1], ac[2][1])
    else
        # swap ac[1,2] with ac[2, 1]
        swap_cols(mat, ac[1][2], ac[2][1])
    end

    return 1

end # function
