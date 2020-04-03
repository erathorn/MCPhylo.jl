function neighbor_joining(dm::Array{Int64,2}, tree::T) where {T<:AbstractNode}
    matrix_size = size(dm)
    n = matrix_size[1]
    leaves = get_leaves(tree)
    # df = DataFrame()
    # leaf_names = [name for leaves.name in leaves]
    # for i in 1:n
    #     df.leaves[i].name = dm[:,i]
    q1 = zeros(Float64, n, n)
    for i = 1:n
        for j = i+1:n
            q1[i, j] =
                q1[j, i] = (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
        end # end for
    end # end for
    index = findmin(q1)[2]
    first_node = leaves[index[1]]
    second_node = leaves[index[2]]
    blub = Node("Blub")
    blub.children = [first_node, second_node]
    blub.nchild = 2
    first_node.mother = blub
    second_node.mother = blub
    first_distance =
        0.5 * dm[index] + (1 / (2 * (n - 2))) * (sum(dm[i, :]) - sum(dm[j, :]))
    first_node.inc_length = first_distance
    second_node.inc_length = dm[index] - first_distance

    dm2 = zeros(Float64, n-1, n-1)
    let j = 1
        for i = 2:n-1
            while  j == index[1] || j == index[2]
                j += 1
            end # end while
            dm2[i,1] = dm2[1,i] =  0.5 * (dm[index[1],j] + dm[index[2],j] - dm[index[1],index[2]])
            j += 1
        end
    end
    for k = 2:n-1
        for l = 2:n-1
            for i = 1:n
                if i == index[1] || i == index[2]
                    continue
                end
                for j = 1:n
                    if j == index[1] || j == index[2]
                        print("I'm here")
                    else
                        dm2[k,l] = dm[i,j]
                        print("update")
                    end
                end
            end
        end
    return dm2
    end

end # end function neighbor_joining
