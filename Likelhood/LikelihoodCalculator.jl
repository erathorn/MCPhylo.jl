module LikelihoodCalculator

using Markdown
include("../Tree/Tree_Module.jl")
using ..Tree_Module
include("../Substitution/SubstitutionMat.jl")
import ..SubstitutionMat: exponentiate_binary

# TODO: Look into Parallelizing it
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node}, pi::Number)#, rates::Vector{Float64})
    root = last(tree_postorder)
    n_c::Int64 = size(root.data)[2]
    for node in tree_postorder
        if node.nchild != 0
            left_daughter = node.child[1]
            right_daughter = node.child[2]
            linc::Float64 = left_daughter.inc_length
            rinc::Float64 = right_daughter.inc_length
            # use the inbounds decorator to enable SIMD
            # SIMD greatly improves speed!!!
            @inbounds for ind in 1:n_c
                #r = rates[ind]
                r = 1.0 # r is going to vary
                left_mat = exponentiate_binary(pi, linc, r)
                right_mat = exponentiate_binary(pi, rinc, r)
                """
                For legacy and debugging purposes:
                #a = left_daughter.data[1,ind]*left_mat[1,1] + left_daughter.data[2,ind]*left_mat[2,1]
                #b = left_daughter.data[2,ind]*left_mat[2,2] + left_daughter.data[1,ind]*left_mat[1,2]
                #c = right_daughter.data[1,ind]*right_mat[1,1] + right_daughter.data[2,ind]*right_mat[2,1]
                #d = right_daughter.data[2,ind]*right_mat[2,2] + right_daughter.data[1,ind]*right_mat[1,2]
                #node.data[1,ind] = a*c
                #node.data[2,ind] = b*d
                """
                a = sum(left_mat.*left_daughter.data[:,ind], dims=1)
                c = sum(right_mat.*right_daughter.data[:,ind], dims=1)
                node.data[:,ind] = a.*c

            end # for
        end # if

    end # for

    # sum the two rows
    return sum(log.(root.data.*[pi, 1.0-pi]))
end # function

end  # module LikelihoodCalculator
