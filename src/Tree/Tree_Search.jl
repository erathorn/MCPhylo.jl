
#################### Safe searches ####################

#"""
# proper Markdown comments are not possible
#
#This part creates functions which enable the search for different nodes in the
#tree. It is possible to look for a node via its name, its binary representation
#or to find the root.
#This functionality can be extended by adding more fields to the nodes and the
#meta programmming part here.

# These functions are safe, in the sense that they properly exit if the node is
# not found.
#"""
for (sym, my_type) in [(:binary, :String), (:name, :String), (:root ,:Bool), (:num, :Int64)]
    # extend the list to look for more fields in the node
    @eval function $(Symbol(string("find_by_$sym")))(tree::T, identifier::$my_type)::T  where T<:AbstractNode
        # create each function and make it so it only accepts the correct type
        local all_nodes = post_order(tree) # make sure all_nodes only belongs to this function
        for node in all_nodes
            if node.$sym == identifier
                # return the node if it is found
                return node
            end # if
        end # for
        # the node is not found. Therefore throw an error!
        throw("The node identified by $identifier is not in the tree.")
    end # function
end


#################### Unsafe searches ####################

"""
    find_num(root::T, num::Int64)  where T<:AbstractNode

Find a node by its number. The function assumes that the node is present in the
tree.

Do not use this function if you are unsure wheter the node is in the tree at all.
"""
function find_num(root::T, num::I)::T  where {T<:AbstractNode, I<:Integer}
    rn = Vector{T}(undef, 1)
    find_num(root, num, rn)
    return rn[1]
end


"""
    find_num(root::T, num::Int64, rn::Vector{T})::Bool  where T<:AbstractNode

Do a post order traversal to find the node corresponding to the `num`.
"""
function find_num(root::T, num::I, rn::Vector{T})::Bool  where {T<:AbstractNode, I<:Integer}
    # if the current node is the correct one store it in rn
    if root.num === num
        rn[1] = root
        found = true
    else
        found = false
    end

    if !found
        # if the node is not yet found continue
        for child in root.children
            found = find_num(child, num,  rn)
        end
    end # if
    return found
end

"""
    find_binary(root::T, bin::String)::T where T<:AbstractNode

Find a node by its binary representation. The function assumes that the node is
present in the tree.

Do not use this function if you are unsure wheter the node is in the tree at all.
"""
function find_binary(root::T, bin::String)::T where T<:AbstractNode
    rv = root
    for i in bin[2:end]
        rv = rv.children[parse(Int64, i)+1]
    end
    rv
end


function find_binary(root::T, bin::String, rn::Vector{T})::Bool  where T<:AbstractNode
    # if the current node is the correct one store it in rn
    println(root.num)
    if root.binary == bin
        rn[1] = root
        found = true
    else
        found = false
    end # if
    if !found
        # if the node is not yet found continue
        for child in root.children
            found = find_binary(child, bin,  rn)
        end # for

    end # if
    return found
end
