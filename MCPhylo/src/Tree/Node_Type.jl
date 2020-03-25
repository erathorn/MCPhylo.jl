
"""
    Node

This data type holds the basic Node structure. The type T is used to specify the type of the data
stored in the node.

* If `nchild` is `0` the Node is a leaf node.
* If `root` is `False` the Node is a child of another node.
* `inc_length` specifies the length of the incomming branch.
* `binary` specifies the path from the root to the Node. `1` and `0` represent left and right turns respectively.
"""

abstract type AbstractNode{T,A,B,I} end

mutable struct Node_cu{T<:Real, A<:AbstractArray,B<:AbstractArray, I<:Integer} <: AbstractNode{T,A,B,I}
    name::String
    data::A#{Float64, 2}
    mother::Union{Node_cu, Missing}
    children::Vector{Node_cu}
    nchild::I
    root::Bool
    scaler::B
    inc_length::T
    binary::String
    num::I
    height::Float64
    IntExtMap::Union{Vector{Int64}, Nothing}
    blv::Union{Vector{Float64}, Nothing}
 end


mutable struct Node{T<:Real, A<:AbstractArray,B<:AbstractArray, I<:Integer} <: AbstractNode{T,A,B,I}
    name::String
    data::A
    mother::Union{Node, Missing}
    children::Vector{Node}
    scaler::B
    nchild::I
    root::Bool
    inc_length::T
    binary::String
    num::I
    height::T
    IntExtMap::Union{Vector{I}, Nothing}
    blv::Union{Vector{T}, Nothing}
    initialized::Bool
end # struct Node

function Node()::Node
        Node{Float64,Array{Float64},Array{Float64},Int64}("no_name", ones(3,3), missing,Vector{Node}(undef, 0) ,ones(3),0,true,0.5,"0",1,0.5,nothing,nothing,false)
end

function Node(name::String; data::A)::Node where A<:AbstractArray
        Node{Float64,A,Array{Float64},Int64}(name, data ,missing, Vector{Node}(undef, 0), ones(3), 0, true, 0.5, "0", 1, 0.5, nothing, nothing, false)
end
