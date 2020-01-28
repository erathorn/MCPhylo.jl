
function RF(tree1::Node, tree2::Node)::Int64

    bt3 = get_bipartitions(tree1)
    bt4 = get_bipartitions(tree2)

    length(bt3)+length(bt4) - 2* length(intersect(bt3, bt4))
end

function get_bipartitions(tree::Node)::Vector{Tuple}
    po_vect= post_order(tree)[1:end-1]
    bt = Vector{Tuple}(undef, length(po_vect))
    Base.Threads.@threads for ind in eachindex(po_vect)
        elem = po_vect[ind]::Node
        inset = Set{Int64}()
        outset = Set{Int64}()
        @simd for elem2 in po_vect
            startswith(elem2.binary, elem.binary) ? push!(inset, elem2.num) : push!(outset, elem2.num)
        end
        @inbounds bt[ind] = (inset, outset)
    end

    bt
end

function BHV_lower(tree1::Node, tree2::Node)
    #res_low = 0.0
    res_upper_1 = 0.0
    res_upper_2 = 0.0
    res_upper_3 = 0.0
    po = post_order(tree1)
    Base.Threads.@threads for node in po
        nom = MCPhylo.find_num(tree2, node.num)
        if !node.root
            if node.mother.num == nom.mother.num
                #res_low += (node.inc_length-nom.inc_length)^2
                res_upper_3 += (node.inc_length-nom.inc_length)^2
            else
                #res_low += nom.inc_length^2
                res_upper_2 += nom.inc_length^2
                #res_low += node.inc_length^2
                res_upper_1 += node.inc_length^2
            end
        end
    end
    res_low = res_upper_1+res_upper_2+res_upper_3
    res_high = (sqrt(res_upper_2)+sqrt(res_upper_1))^2+res_upper_3
    sqrt(res_low), sqrt(res_high)
end
