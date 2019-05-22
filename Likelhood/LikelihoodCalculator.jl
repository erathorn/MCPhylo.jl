module LikelihoodCalculator

using Markdown
include("../Tree/Tree_Module.jl")
using ..Tree_Module
include("../Substitution/SubstitutionMat.jl")
using ..SubstitutionMat: exponentiate_binary

<<<<<<< HEAD
# TODO: Look into Parallelizing it
=======
>>>>>>> 394b950e973d95ed7cd59984eed0aacda1aab011
"""
    FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})

This function calculates the log-likelihood of an evolutiuonary model using the
Felsensteins pruning algorithm.
"""
function FelsensteinFunction(tree_postorder::Vector{Node}, pi::Float64, rates::Vector{Float64})
    for node in tree_postorder
        if node.nchild != 0
            left_daughter = node.child[1]
            right_daugher = node.child[2]
            for ind in 1:size(left_daughter.data.size[2])
                r = rates[ind]
                left_mat = SubstitutionMat.exponentiate_binary(pi, left_daughter.inc_length, r)
                right_mat = SubstitutionMat.exponentiate_binary(pi, right_daughter.inc_length, r)
                a = left_daughter.data[ind][1]*left_mat[1][1] + left_daughter.data[ind][2]*left_mat[2][1]
                b = left_daughter.data[ind][2]*left_mat[2][2] + left_daughter.data[ind][1]*left_mat[1][2]
                c = right_daughter.data[ind][1]*right_mat[1][1] + right_daughter.data[ind][2]*right_mat[2][1]
                d = right_daughter.data[ind][2]*right_mat[2][2] + right_daughter.data[ind][1]*right_mat[1][2]
                node.data[ind][1] = a+c
                node.data[ind][2] = b+d
            end # for
        end # if
    end # for
    root = last(tree_postorder)
    # sum the two rows
    return sum(log.(sum(root.data.*[pi, 1.0-pi], dims=2)))
end # function

end  # module LikelihoodCalculator
