include("./MCPhylo/src/MCPhylo.jl")

"""
neighbor_joining(dm::ArrayFloat64,2}, tree::T)

This function returns a phylogenetic tree by using neighbor joining
based on a given distance matrix
"""
# Tree structure, node naming, node numbers needed?, example tree?, let blocks?
# use RF function to check topology
function neighbor_joining(dm::Array{Float64,2}, leaves::Vector{T}) where {T>:AbstractNode}

    n = size(dm)[1]
    while n > 2
        q1 = zeros(Float64, n, n)
        for i = 1:n
            for j = i+1:n
                q1[i, j] =
                    q1[j, i] =
                        (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
            end # end for
        end # end for

        println(q1)

        count = 1
        index = findmin(q1)[2]
        first_node = leaves[index[2]]
        second_node = leaves[index[1]]
        new_node = Node()
        new_node.name = "Node_$count"
        new_node.children = [first_node, second_node]
        new_node.nchild = 2

        println(new_node)

        first_node.mother = new_node
        second_node.mother = new_node
        if index[1] < index[2]
            deleteat!(leaves, [index[1], index[2]])
        else
            deleteat!(leaves, [index[2], index[1]])
        end

        insert!(leaves, 1, new_node)

        dummy = 0.5 * dm[index] + (1 / (2 * (n - 2))) * (sum(dm[index[2], :]) - sum(dm[index[1], :]))

        println(dummy)

        first_node.inc_length = dummy
        second_node.inc_length = dm[index] - dummy

        println(first_node)
        println(second_node)

        next_dm = zeros(Float64, n - 1, n - 1)
        let j = 1
            for i = 2:n-1
                while j == index[1] || j == index[2]
                    j += 1
                end # end while
                next_dm[i, 1] =
                    next_dm[1, i] =
                        0.5 * (
                            dm[index[1], j] + dm[index[2], j] -
                            dm[index[1], index[2]]
                        )
                j += 1
            end # end for
        end # end let

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

        let z = 1
            for k = 2:n-1
                for l = 2:n-1
                    next_dm[k, l] = entries[z]
                    z += 1
                end   # end for
            end # end for
        end # end let
        n -= 1 # update matrix dimension
        dm = next_dm
        #println(dm)
        count += 1
    end # end while
    root = last(leaves)
    root.children = [leaves[1]]
    root.nchild = 1
    leaves[1].mother = root

    return root
end # end function neighbor_joining


function neighbor_joining(dm::Array{Float64,2})

    n = size(dm)[1]
    leaves = AbstractNode[]
    for i = 1:n
        append!(leaves, Node("leaf_$i"))
    end
    while n > 2
        q1 = zeros(Float64, n, n)
        for i = 1:n
            for j = i+1:n
                q1[i, j] =
                    q1[j, i] =
                        (n - 2) * dm[i, j] - sum(dm[i, :]) - sum(dm[j, :])
            end # end for
        end # end for
        count = 1
        index = findmin(q1)[2]

        first_node = leaves[index[1]]
        second_node = leaves[index[2]]
        new_node = Node()
        new_node.name = "Node_$count"
        new_node.children = [first_node, second_node]
        new_node.nchild = 2
        first_node.mother = new_node
        second_node.mother = new_node
        deleteat!(leaves, [index[1], index[2]])
        insert!(leaves, 1, new_node)

        first_node.inc_length =
            0.5 * dm[index] +
            (1 / (2 * (n - 2))) * (sum(dm[index[1], :]) - sum(dm[index[2], :]))
        second_node.inc_length = dm[index] - first_node.inc_length

        next_dm = zeros(Float64, n - 1, n - 1)
        let j = 1
            for i = 2:n-1
                while j == index[1] || j == index[2]
                    j += 1
                end # end while
                next_dm[i, 1] =
                    next_dm[1, i] =
                        0.5 * (
                            dm[index[1], j] + dm[index[2], j] -
                            dm[index[1], index[2]]
                        )
                j += 1
            end # end for
        end # end let

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

        let z = 1
            for k = 2:n-1
                for l = 2:n-1
                    next_dm[k, l] = entries[z]
                    z += 1
                end   # end for
            end # end for
        end # end let
        n -= 1 # update matrix dimension
        dm = next_dm
        count += 1
    end # end while
    root = last(leaves)
    root.children = [leaves[1]]
    root.nchild = 1
    leaves[1].mother = root

    return root
end # end function neighbor_joining

function test_nj()

    dm = [0.0 5.0 9.0 9.0 8.0;5.0 0.0 10.0 10.0 9.0;9.0 10.0 0.0 8.0 7.0;
    9.0 10.0 8.0 0.0 3.0;8.0 9.0 7.0 3.0 0.0]
    a = Node()
    a.name = "a"
    b = Node()
    b.name = "b"
    c = Node()
    c.name = "c"
    d = Node()
    d.name = "d"
    e = Node()
    e.name = "e"
    leaves = [a, b, c, d, e]
    neighbor_joining(dm, leaves)

end # end function test_nj
