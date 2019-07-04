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
function FelsensteinFunction(tree_postorder::Vector{Node}, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Float64

    for node in tree_postorder
        if node.nchild !== 0
            CondLikeInternal(node, pi_, rates, n_c)
        end # if
    end # for

    # sum the two rows
    rdata::Array{Float64,2}=last(tree_postorder).data
    res::Float64 = 0.0
    _pi_::Float64 = log(1.0-pi_)
    _lpi_::Float64 = log(pi_)
    @inbounds for ind in 1:n_c
        res +=(log(rdata[1,ind])+ _lpi_) + (log(rdata[2,ind])+ _pi_)

        #rdata[1, ind] *pi_
        #rdata[2, ind] *=_pi_
    end

    return res#sum(log.(rdata.*[pi_, 1.0-pi_]))
end # function

function CondLikeInternal(node::Node, pi_::Number, rates::Vector{Float64}, n_c::Int64)::Nothing
    @assert size(node.child)[1] == 2
    @assert size(rates)[1] == n_c
    left_daughter::Node = node.child[1]
    right_daughter::Node = node.child[2]
    linc::Float64 = left_daughter.inc_length
    rinc::Float64 = right_daughter.inc_length
    left_daughter_data::Array{Float64,2} = left_daughter.data
    right_daughter_data::Array{Float64,2} = right_daughter.data

    # use the inbounds decorator to enable SIMD
    # SIMD greatly improves speed!!!
    @inbounds for ind in 1:n_c
        r::Float64 = rates[ind]
        left_mat::Array{Float64,2} = exponentiate_binary(pi_, linc, r)
        right_mat::Array{Float64,2} = exponentiate_binary(pi_, rinc, r)
        #lv1::Float64 = left_daughter_data[1,ind]
        #lv2::Float64 = left_daughter_data[2,ind]
        #rv1::Float64 = right_daughter_data[1,ind]
        #rv2::Float64 = right_daughter_data[2,ind]

        a::Float64 = left_daughter_data[1,ind]*left_mat[1,1] + left_daughter_data[2,ind]*left_mat[2,1]
        b::Float64 = left_daughter_data[1,ind]*left_mat[1,2] + left_daughter_data[2,ind]*left_mat[2,2]
        c::Float64 = right_daughter_data[1,ind]*right_mat[1,1] + right_daughter_data[2,ind]*right_mat[2,1]
        d::Float64 = right_daughter_data[1,ind]*right_mat[1,2] + right_daughter_data[2,ind]*right_mat[2,2]

        node.data[1,ind] = a*c
        node.data[2,ind] = b*d



    end # for
end # function


end  # module LikelihoodCalculator
