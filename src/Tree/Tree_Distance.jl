# TODO: check trees are on the same leave set

"""
    RF(tree1::T, tree2::T)::Int64 where T <:GeneralNode

Calculate the Robinson-Foulds distance between the two trees.
In its current form the function assumes the trees have identical leave sets.
"""
function RF(tree1::T, tree2::T)::Int64 where T <: GeneralNode
    bt3 = get_bipartitions(tree1)
    bt4 = get_bipartitions(tree2)
    length(bt3)+length(bt4) - 2* length(intersect(bt3, bt4))
end

"""
    get_bipartitions(tree::T)::Vector{Tuple} where T <:GeneralNode

Get a vector of all bipartions of `tree`. The resulting vector contains Tuples
of sets representing the bipartions.
"""
function get_bipartitions(tree::T)::Vector{Tuple} where T <:GeneralNode
    po_vect= post_order(tree)[1:end-1]
    bt = Vector{Tuple}(undef, length(po_vect))
    Base.Threads.@threads for ind in eachindex(po_vect)
        elem = po_vect[ind]::T
        inset = Set{Int64}()
        outset = Set{Int64}()
        @simd for elem2 in po_vect
            startswith(elem2.binary, elem.binary) ? push!(inset, elem2.num) : push!(outset, elem2.num)
        end # for
        @inbounds bt[ind] = (inset, outset)
    end # for
    bt
end

"""
    BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:GeneralNode

This function calculates the lower and upper bounds of the geodesic in the
Billera-Holmes-Vogtman space.
"""
function BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:GeneralNode
    # res_upper_1 = 0.0
    # res_upper_2 = 0.0
    # res_upper_3 = 0.0
    # po = post_order(tree1)
    # Base.Threads.@threads for node in po
    #     nom = MCPhylo.find_num(tree2, node.num)
    #     if !node.root
    #         if node.mother.num == nom.mother.num
    #             #res_low += (node.inc_length-nom.inc_length)^2
    #             res_upper_3 += (node.inc_length-nom.inc_length)^2
    #         else
    #             #res_low += nom.inc_length^2
    #             res_upper_2 += nom.inc_length^2
    #             #res_low += node.inc_length^2
    #             res_upper_1 += node.inc_length^2
    #         end
    #     end
    # end
    # res_low = res_upper_1+res_upper_2+res_upper_3
    # res_high = (sqrt(res_upper_2)+sqrt(res_upper_1))^2+res_upper_3
    # sqrt(res_low), sqrt(res_high)
    #throw("not here")
    res_upper_1 = Threads.Atomic{Float64}(0.0)
    res_upper_2 = Threads.Atomic{Float64}(0.0)
    res_upper_3 = Threads.Atomic{Float64}(0.0)
    po::Vector{T} = post_order(tree1)
    Base.Threads.@threads for node in po[1:end-1]
        nom = find_num(tree2, node.num)
        if node.mother.num == nom.mother.num
            Threads.atomic_add!(res_upper_3, (node.inc_length-nom.inc_length)^2)
            #res_upper_3 += (node.inc_length-nom.inc_length)^2
        else
            Threads.atomic_add!(res_upper_2 , nom.inc_length^2)
            #res_upper_2 += nom.inc_length^2
            Threads.atomic_add!(res_upper_1, node.inc_length^2)
            #res_upper_1 += node.inc_length^2
        end # if

    end # for
    res_low::Float64 = res_upper_1[]+res_upper_2[]+res_upper_3[]
    #res_low::Float64 = res_upper_1+res_upper_2+res_upper_3
    res_high::Float64 = (sqrt(res_upper_2[])+sqrt(res_upper_1[]))^2+res_upper_3[]
    #res_high::Float64 = (sqrt(res_upper_2)+sqrt(res_upper_1))^2+res_upper_3
    sqrt(res_low), sqrt(res_high)
end
