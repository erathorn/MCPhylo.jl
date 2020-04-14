include("./MCPhylo/src/MCPhylo.jl")

"""
neighbor_joining(dm::Array{Int64,2}, Array{String,1})

This function returns a phylogenetic tree by using neighbor joining
based on a given distance matrix and an array of leaf names
"""
function neighbor_joining(dm::Array{Float64,2}, leaf_names::Array{String,1})

    n = size(dm)[1]
    leaves = AbstractNode[]
    for leaf in leaf_names
        new_leaf = Node()
        new_leaf.name = leaf
        push!(leaves, new_leaf)
    end # end for
    # count for node names
    count = 1
    while n > 2
        # build Q1 matrix
        q1 = zeros(Float64, n, n)
        for i = 1:n
            for j = i+1:n
                q1[i, j] =
                    q1[j, i] =
                        (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
            end # end for
        end # end for
        # locate minimum of Q1
        index = findmin(q1)[2]

        first_node = leaves[index[2]]
        second_node = leaves[index[1]]
        new_node = Node()
        new_node.name = "Node_$count"
        count += 1
        new_node.children = [first_node, second_node]
        new_node.nchild = 2
        first_node.mother = new_node
        second_node.mother = new_node

        println(new_node.name)
        println(first_node.name)
        println(second_node.name)
        # update array with leaf names
        deleteat!(leaves, [index[2]])
        deleteat!(leaves, [index[1] - 1])
        insert!(leaves, 1, new_node)
        # calculate distance
        first_node.inc_length =
            0.5 * dm[index] +
            (1 / (2 * (n - 2))) * (sum(dm[index[2], :]) - sum(dm[index[1], :]))
        second_node.inc_length = dm[index] - first_node.inc_length

        println([first_node.inc_length, second_node.inc_length])

        # initalize next distance matrix
        next_dm = zeros(Float64, n - 1, n - 1)
        # fill first row and column of next distance matrix
        j = 1
        for i = 2:n-1
            while j == index[1] || j == index[2]
                j += 1
            end # end while
            next_dm[i, 1] =
                next_dm[1, i] =
                    0.5 *
                    (dm[index[1], j] + dm[index[2], j] - dm[index[1], index[2]])
            j += 1
        end # end for
        # copy values of last distance matrix to fill the next one
        entries = Float64[]
        for i = 1:n
            if i == index[1] || i == index[2]
                continue
            end # end if
            for j = 1:n
                if j == index[1] || j == index[2]
                    continue
                else
                    append!(entries, dm[i, j])
                end # end else
            end # end for
        end # end for
        # fill remaining cells of next distance matrix
        z = 1
        for k = 2:n-1
            for l = 2:n-1
                next_dm[k, l] = entries[z]
                z += 1
            end   # end for
        end # end for
        n -= 1 # update matrix dimension
        dm = next_dm
        # add final (root) node to tree
        if n == 2
            root = last(leaves)
            root.root = true
            child_node = leaves[1]
            root.children = [child_node]
            root.nchild = 1
            child_node.mother = root
            child_node.inc_length = dm[index]
            println(child_node.inc_length)
            return root
        end # end if
    end # end while
end # end function neighbor_joining

"""
neighbor_joining(dm::Array{Int64,2}, Array{String,1})

This function returns a phylogenetic tree by using neighbor joining
based on a given distance matrix
"""
function neighbor_joining(dm::Array{Float64,2})

    n = size(dm)[1]
    leaves = AbstractNode[]
    # build array of dummy leaves
    for i = 1:n
        new_leaf = Node()
        new_leaf.name = "leaf_$i"
        push!(leaves, new_leaf)
    end # end for
    # count for node names
    count = 1
    while n > 2
        # build Q1 matrix
        q1 = zeros(Float64, n, n)
        for i = 1:n
            for j = i+1:n
                q1[i, j] =
                    q1[j, i] =
                        (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
            end # end for
        end # end for
        # locate minimum of Q1
        index = findmin(q1)[2]

        first_node = leaves[index[2]]
        second_node = leaves[index[1]]
        new_node = Node()
        new_node.name = "Node_$count"
        count += 1
        new_node.children = [first_node, second_node]
        new_node.nchild = 2
        first_node.mother = new_node
        second_node.mother = new_node

        println(new_node.name)
        println(first_node.name)
        println(second_node.name)
        # update array with leaves
        deleteat!(leaves, [index[2]])
        deleteat!(leaves, [index[1] - 1])
        insert!(leaves, 1, new_node)
        # calculate distance
        first_node.inc_length =
            0.5 * dm[index] +
            (1 / (2 * (n - 2))) * (sum(dm[index[2], :]) - sum(dm[index[1], :]))
        second_node.inc_length = dm[index] - first_node.inc_length

        println([first_node.inc_length, second_node.inc_length])

        # initalize next distance matrix
        next_dm = zeros(Float64, n - 1, n - 1)
        # fill first row and column of next distance matrix
        j = 1
        for i = 2:n-1
            while j == index[1] || j == index[2]
                j += 1
            end # end while
            next_dm[i, 1] =
                next_dm[1, i] =
                    0.5 *
                    (dm[index[1], j] + dm[index[2], j] - dm[index[1], index[2]])
            j += 1
        end # end for

        # copy values of last distance matrix to fill the next one
        entries = Float64[]
        for i = 1:n
            if i == index[1] || i == index[2]
                continue
            end # end if
            for j = 1:n
                if j == index[1] || j == index[2]
                    continue
                else
                    append!(entries, dm[i, j])
                end # end else
            end # end for
        end # end for
        # fill remaining cells of next distance matrix
        z = 1
        for k = 2:n-1
            for l = 2:n-1
                next_dm[k, l] = entries[z]
                z += 1
            end   # end for
        end # end for
        n -= 1 # update matrix dimension
        dm = next_dm
        # add final (root) node to tre
        if n == 2
            root = last(leaves)
            root.root = true
            child_node = leaves[1]
            root.children = [child_node]
            root.nchild = 1
            child_node.mother = root
            child_node.inc_length = dm[index]
            println(child_node.inc_length)
            return root
        end # end if
    end # end while
end # end function neighbor_joining
