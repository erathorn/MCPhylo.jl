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
    for node in tree_postorder
        if node.nchild != 0
            left_daughter = node.child[1]
            right_daughter = node.child[2]
            for ind in 1:size(left_daughter.data)[2]
                #r = rates[ind]
                r = 1.0
                left_mat = exponentiate_binary(pi, left_daughter.inc_length, r)
                right_mat = exponentiate_binary(pi, right_daughter.inc_length, r)

                a = (left_daughter.data[1,ind]*left_mat[1,1]) + (left_daughter.data[2,ind]*left_mat[1,2])
                b = left_daughter.data[2,ind]*left_mat[2,2] + left_daughter.data[1,ind]*left_mat[2,1]
                c = right_daughter.data[1,ind]*right_mat[1,1] + right_daughter.data[2,ind]*right_mat[1,2]
                d = right_daughter.data[2,ind]*right_mat[2,2] + right_daughter.data[1,ind]*right_mat[2,1]
                node.data[1,ind] = a*c
                node.data[2,ind] = b*d
            end # for
        end # if

    end # for
    root = last(tree_postorder)
    # sum the two rows
    return log(sum(exp.(sum(log.(root.data.*[pi, 1.0-pi])))))
end # function

end  # module LikelihoodCalculator
