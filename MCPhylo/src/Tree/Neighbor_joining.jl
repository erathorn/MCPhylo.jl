function neighbor_joining(dm::Array{Int64, 2}, tree::T) where T <:AbstractNode
    matrix_size = size(dm)
    n = matrix_size[1]
    leaves = get_leaves(tree)
    df = DataFrame()
    leaf_names =
    for i in 1:n
        df.leaves[i].name = dm[:,i]
    q1 = zeros(Float64, n, n)
    for i in 1:n
        for j in i+1:n
            q1[i,j] = (n-2) * dm[i,j] - sum(dm[i,:]) - sum(dm[j,:])
            q1[j,i] = q1[i,j]
        end # end for
    end # end for
    index = findmin(q1)[2]

    first_node = leaves[index[1]]
    second_node = leaves[index[2]]


    df = DataFrame(Names = dm[:,1], Age = dm[:,2], Gender = dm[:,3], Job = dm[:,4], Nationality = dm[:,5])
end # end function neighbor_joining
