
"""
    Node

This data type holds the basic Node structure. The type T is used to specify the type of the data
stored in the node.

* If `nchild` is `0` the Node is a leaf node.
* If `root` is `False` the Node is a child of another node.
* `inc_length` specifies the length of the incomming branch.
* `binary` specifies the path from the root to the Node. `1` and `0` represent left and right turns respectively.
"""

abstract type AbstractNode end

mutable struct Node_cu{T<:Real, A<:AbstractArray,B<:AbstractArray, I<:Integer} <: AbstractNode
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


mutable struct Node{T<:Real, A<:Real,B<:Real, I<:Integer} <: AbstractNode
    name::String
    data::Array{A,2}
    mother::Union{Node{T,A,B,I}, Missing}
    children::Vector{Node{T,A,B,I}}
    scaler::Array{B,2}
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


#const NodeS = Node{Float64,Array{Float64,2},Array{Float64,1},Int64}

function Node()::Node
        Node{Float64,Float64,Float64,Int64}("no_name", ones(3,3), missing,Vector{Node{Float64,Float64,Float64,Int64}}(undef, 0) ,ones(1,3),0,true,0.5,"0",1,0.5,nothing,nothing,false)
end

function Node(name::String; data::Array{A,2})::Node where A<:Real
        Node{Float64,A,Float64,Int64}(name, data ,missing, Vector{Node{Float64,A,Float64,Int64}}(undef, 0), ones(3,1), 0, true, 0.5, "0", 1, 0.5, nothing, nothing, false)
end



#################### Base functionality ####################

Base.:(==)(x::T, y::T) where T<:AbstractNode = x.num == y.num
Base.size(x::T) where T<:AbstractNode = size(post_order(x))
Base.length(x::T) where T<:AbstractNode = x.nchild

function Base.summary(io::IO, d::N) where N <: AbstractNode
    summary(io, d.name)
end


function Base.show(io::IO, d::N) where N <: AbstractNode
    print(io, "Tree with root:\n")
    show(io, d.name)
    if d.initialized
        print(io, "\nLength:\n")
        show(io, "text/plain", tree_length(d))
        print(io, "\nHeight:\n")
        show(io, "text/plain", tree_height(d))
        print(io, "\nNumber of leave nodes:\n")
        show(io, "text/plain",length(get_leaves(d)))
    else
        print(io, "\nLength:\n")
        show(io, "text/plain", 0)
        print(io, "\nHeight:\n")
        show(io, "text/plain", 0)
        print(io, "\nNumber of leave nodes:\n")
        show(io, "text/plain",0)
    end
end

function showall(io::IO, d::N) where N <: AbstractNode
  show(io, d)
  print(io, "\nNode:\n")
  show(io, "text/plain", d.name)
  print(io, "\n#children:\n")
  show(io, d.nchild)
  print(io, "\nbinary:\n")
  show(io, d.binary)
end
