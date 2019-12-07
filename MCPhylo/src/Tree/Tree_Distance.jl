
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
