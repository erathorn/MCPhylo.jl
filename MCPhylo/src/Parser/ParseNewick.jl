
# that's very far from ideal, but atom and I don't understand each other otherwise
include("../MCPhylo.jl")

"""
    ParseNewick(filename::String)

This function parses a Newick file
"""
function ParseNewick(filename::String)
    content = load_newick(filename)
    if !is_valid_newick_string(content)
        throw("$filename is not a Newick file!")
    end # if
    content = strip(content)

    #the parse thing
    # TODO: actually useful part of the code goes here


end


"""
    load_newick(filename::String)

This function loads a newick from file
"""
# TODO: rn it assumed that there are no extra line breaks (\n\n) and there is only one tree pro file. That should be fixed.

function load_newick(filename::String)
    open(filename, "r") do file
        global content = readlines(file)
    end
    content[1]
end

"""
    is_valid_newick_string(newick::String)

This function checks if the given string is valid: is the brackets number matches and if the string ends with ";"
"""

function is_valid_newick_string(newick::String)
    # Step one: does the stripped string ends with ';'
    if endswith(strip(newick),";")
        # Step two: check for the equal amount of brackets
        bracket_level = 0
        for char in newick
            if char == '('
                bracket_level += 1
            elseif char == ')'
                    bracket_level -= 1
            end # elseif
        end # for
        if bracket_level != 0
            return false
        end # if
    else # same level as the endswith statement
        return false
    end # else
    return true
end

# the grand line between functions which are we more certain in and experimental functions on which we are still working on

###########

function parsing_the_newick(newick::String,current_node::Any, cur_loc::Int)

    if current_node == nothing
        cur_loc = 1
        count = 0
        current_node = Node()
        println("Starting up!")
    end # setting things up

    if newick[cur_loc] == '('

        cur_loc+=1
        println("Current location is updated, it's ", cur_loc)

        if newick[cur_loc] == '('
            # YOUR RECURSION IS HERE

            left_child = parsing_the_newick(string(SubString(newick,cur_loc)),current_node,cur_loc)
            println("Right now we are working with the string ", string(SubString(newick,cur_loc)))
            cur_loc+=1
            add_child!(current_node,left_child)
        end # if (the recursive call one)
        node_boarder = match(r"[();,]",string(SubString(newick,cur_loc))).offset
        println("The node boarder is ",node_boarder)
        name,length = parse_name_length(string(SubString(newick,cur_loc,node_boarder+1)))
        if name != "no_name"
            println("This node has a name and length")
            left_child = Node()
            left_child.name = name
            left_child.inc_length = length
            left_child.mother = current_node
            current_node.nchild +=1
            left_child.num = count
            count++
            cur_loc = node_boarder
            print(left_child.mother)
        else
            println("You're nameless/And that's your future from now on.")
        end # if_else
    end # first if
    return current_node
end # the function


"""
    parse_name_length(newick::String)

This function parses two optional elements of the tree, name and length. In case, when neither of this is provided, empty string and nothing are return
"""

function parse_name_length(newick::String)
    if occursin(':',newick)
        name, length = split(newick,':')
        return string(name), parse(Float64, length)
    end # if
    "no_name", nothing
end




parsing_the_newick("((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;",nothing,1)





       # // try to find a valid $name next at the current cursor position //
       # if $name = $newick_string [ $cursor ] . regex match ("^[0-9A-Za-z_|]+") {
       #     // Then the $left node is a leaf node: parse the leaf data //
       #     $this.left = new Object;
       #     $this.left.name = $name;
       #     $this.left.serial = $count;
       #     $count ++;
       #     $this.left.parent = $this;
       #     $cursor += length of $name;
       #     // move cursor to position after matched name.
       # }







# TODO: rewrite this one
#
# function make_node(newick::String)
#     parts = split(newick, ')')
#     if length(parts) == 1
#         label = newick
#         children = Vector{Node}(undef, 0)
#         name, inc_length = parse_name_length(label)
#         return Node{Float64,Array{Float64,2},Array{Float64},Int64}(name,ones(3,3),missing,children,ones(3),3,false,inc_length,"0",1,0.5,nothing,nothing,true)
#     else
#         # TODO: why can't we use length? check
#         len_minus_one = size(parts,1)-1
#         children = list(parse_siblings(join(parts[len_minus_one],')')[2:size(parts,1)]))
#         label = parts[len_minus_one]
#     end #if
#     name, inc_length = parse_name_length(label)
#     parent = Node(name,ones(3,3),missing,children,length(children),false,inc_length,"0",1,0.5,nothing,nothing,true)
#     for x in children
#         x.mother=parent
#     end #for
#     return parent
# end #function
#


function parse_siblings(newick::String)
    bracket_lvl = 0
    current = []
    for c in  (newick * ',')
        if c == ','
            if bracket_lvl == 0
                yield(make_node(join(current,"")))
                current = []
            else
                if c == '('
                    bracket_lvl += 1
                elseif c == ')'
                    bracket_lvl -= 1
                end #elseif
            end #if
            push!(current,c)
        end #if
    end #for
end #function
