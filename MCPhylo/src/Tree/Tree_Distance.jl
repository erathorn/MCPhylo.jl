
function RF(tree1::Node, tree2::Node)
    bt1 = get_bipartitions(tree1)
    bt2 = get_bipartitions(tree2)
    length(bt1)+length(bt2) - 2* length(intersect(bt1, bt2))
end

function get_bipartitions(tree::Node)::Vector
    bt = []
    binlist::Array{String} = []
    numlist::Array{Int64} = []
    for node in post_order(tree)[1:end-1]
        push!(binlist, node.binary)
        push!(numlist, node.num)
    end
    for (ind, elem) in enumerate(binlist)
        inset = Set{Int64}()
        outset = Set{Int64}()
        for (ind2, elem2) in enumerate(binlist)
            if startswith(elem2, elem)
                push!(inset, numlist[ind2])
            else
                push!(outset, numlist[ind2])
            end
        end
        push!(bt, (inset, outset))
    end
    bt
end
