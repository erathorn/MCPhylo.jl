
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

mutable struct GeneralNode{S<: AbstractString, R<:Real, A<:AbstractArray{<:Real},
                    C<:AbstractArray{<:Real}, I<:Integer, T<: AbstractString, B<:Bool} <: AbstractNode
    name::S
    data::A
    mother::Union{GeneralNode{S,R,A,C,I,T,B}, Missing}
    children::Vector{GeneralNode{S,R,A,C,I,T,B}}
    scaler::C
    nchild::I
    root::B
    inc_length::R
    binary::T
    num::I
    height::R
    IntExtMap::Vector{I}
    blv::Vector{R}
    initialized::B
    stats::Dict{String, Float64}
end # struct Node

const Node = GeneralNode{String, Float64, Array{Float64, 2}, Array{Float64, 2}, Int64, String, Bool}
const Node_cu = GeneralNode{String, Float64, CuArray{Float64}, CuArray{Float64}, Int64, String, Bool}


function Node()::Node
        Node("no_name", ones(3,3), missing,Node[] ,ones(1,3),0,true,1.0,"0",1,1.0,Int64[],Float64[],false,  Dict{String, Float64}())
end


function Node(name::String; data::Array{A,2}=ones(2,3))::Node where A<:Real
        Node(name, data ,missing, Node[], ones(3,2), 0, true, 1.0, "0", 1, 1.0, Int64[], Float64[], false,  Dict{String, Float64}())
end


#################### Base functionality ####################

Base.:(==)(x::T, y::T) where T<:GeneralNode = x.num == y.num
Base.size(x::T) where T<:GeneralNode = size(post_order(x))
Base.length(x::T) where T<:GeneralNode = x.nchild

function Base.summary(io::IO, d::N) where N <: GeneralNode
    summary(io, d.name)
end


function Base.show(io::IO, d::N) where N <: GeneralNode
    print(io, "Tree with root:\n")
    show(io, d.name)
    if d.initialized
        print(io, "\nLength:\n")
        show(io, "text/plain", tree_length(d))
        print(io, "\nHeight:\n")
        show(io, "text/plain", tree_height(d))
    else
        print(io, "\nLength:\n")
        show(io, "text/plain", 0)
        print(io, "\nHeight:\n")
        show(io, "text/plain", 0)
        print(io, "\nNumber of leave nodes:\n")
        show(io, "text/plain",0)
    end
end

function showall(io::IO, d::N) where N <: GeneralNode
  show(io, d)
  print(io, "\nNode:\n")
  show(io, "text/plain", d.name)
  print(io, "\n#children:\n")
  show(io, d.nchild)
  print(io, "\nbinary:\n")
  show(io, d.binary)
end
