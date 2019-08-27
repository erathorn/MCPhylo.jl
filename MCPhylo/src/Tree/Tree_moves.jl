
#TODO: RSPR
"""
    NNI!(root::Node)

This function performs an inplace nearest neighbour interchange operation on the
tree which is supplied.
"""
function NNI!(root::Node)

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
    NNI(mat::Array{Float64,2})::nothing

documentation
"""
function NNI!(root::Node, target::Node)::Int64

    if target.nchild == 0 || target.root
        return 0
        #end
    end

    mother::Node = get_mother(root, target)

    if rand([1,2]) == 1
        child1 = remove_child!(mother, true)
        child2 = remove_child!(mother, false)

        gchild1 = remove_child!(target, true)
        gchild2 = remove_child!(target, false)
    else
        child1 = remove_child!(mother, false)
        child2 = remove_child!(mother, true)
        gchild1 = remove_child!(target, true)
        gchild2 = remove_child!(target, false)
    end # if

    add_child!(target, child1, true)
    add_child!(target, gchild1, false)
    add_child!(mother, child2, true)
    add_child!(mother, gchild2, false)

    set_binary!(root)


    return 1

end # function



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

    if rand([true,false])
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
