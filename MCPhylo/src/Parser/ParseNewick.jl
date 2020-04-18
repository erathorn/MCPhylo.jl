
# include("../Tree/Tree_Traversal.jl")
# #include("../MCPhylo.jl")
#include("../Tree/Node_Type.jl")
#include("../Tree/Tree_Basics.jl")



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

#DISCLAIMER: this function is heavily inspired by the pseudocode provided by https://eddiema.ca/2010/06/25/parsing-a-newick-tree-with-recursive-descent/

#currently we shrink the string instead of following it with the cursor. Cursor should be less time complex, but string is more visual
# TODO: rewrite to the cursor version, remove all exessive println'es, etc

#the possible alternative parsing method, who knows
function testing_new_strat(newick::String, current_node::Any, count::Integer)

    if current_node == nothing
        count = 0
        current_node = Node()
    end #if
    println("it begins")
    println("current newick: ",newick)

    while true
         # if newick[1] == '('
         #     #this deals with the open parenthesis
         #     newick = SubString(newick,2)

            if newick[1] == '(' #TODO
                #this is recursion; if we're looking at an internal node this should happen
                println("ALERT: recursion needed")
                childs_section = match(r"\(([^()]|(?R))*\)",newick) #should return only the descendants of the current child node, check https://regex101.com/r/lF0fI1/1 for proof of this regex working the way it should
                childs_section = childs_section.match
                index=findlast(")",childs_section)[1]+1
                the_rest = SubString(childs_section,index)
                newick = the_rest
                childs_section = SubString(childs_section,2,findlast(')',childs_section)-1)
                cur_child = Node()
                add_child!(current_node,testing_new_strat(string(childs_section),cur_child,count)) #where the magic happens
                if newick[1] == ':'
                    node_boarder = match(r"[();,]",newick).offset
                    name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
                    current_node.inc_length = length
                    current_node.num = count
                    newick = SubString(newick,node_boarder)
                end #if
            end #if

            if occursin(r"^[0-9A-Za-z_|]+",string(newick[1])) || newick[1] == ':'
                #this should happen if the current node's a leaf node
                node_boarder = match(r"[();,]",newick).offset
                name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
                cur_child = Node()
                cur_child.name = name
                cur_child.inc_length = length
                cur_child.num = count
                add_child!(current_node,cur_child)
                newick = SubString(newick,node_boarder)
            end #if
            if newick[1] == ',' #we don't need these i don't think, just gotta look at the next thing
                newick = SubString(newick,2)
            end #if
            if newick[1] == ""
                #the third possibility; should just return the current node, move out of recursion, etc
                return current_node
            end #if
        end #while
    end #function

function parsing_the_newick(newick::String,current_node::Any,count::Integer)

    if current_node == nothing
        count = 0
        current_node = Node()
        println()
        println()
        println()
    end # setting things up

    # the child case

    if newick[1] == '(' || newick[1] == ','

        newick = SubString(newick,2)

        if newick[1] == '('
            # YOUR RECURSION IS HERE
            child = parsing_the_newick(string(newick),current_node,count)
            add_child!(current_node,child)
        end # if (the recursive call one)

        node_boarder = match(r"[();,]",newick).offset
        println("Trying to detect the name and length from the ", string(SubString(newick,1,node_boarder-1)))
        name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))

        child = Node()
        child.name = name
        child.inc_length = length
        child.num = count
        count+=1
        add_child!(current_node,child)
        newick = SubString(newick,node_boarder)
end # if "("

if newick[1] == ')'
        newick = string(SubString(newick,2))
        println("Continiue to parse... ", newick)
end # if with the ")" bracket

if occursin(r"^[0-9A-Za-z_|]+",string(newick[1])) || newick[1] == ':'
        println("We're getting the information of the current node")
        node_boarder = match(r"[();,]",newick).offset
        println("Trying to detect the name and length from the ", string(SubString(newick,1,node_boarder-1)))
        name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
        current_node.name = name
        current_node.inc_length = length
        current_node.num = count
        count+=1
        newick = SubString(newick,node_boarder)
        println("Information about the current node was written down. Continue to parse... ",newick)
end # if length one
    if newick[1] == ';'
        println("We have reached the end!")
    end # the last one
return newick, current_node
end # the function



"""
    parse_name_length(newick::String)

This function parses two optional elements of the tree, name and length. In case, when neither of this is provided, empty string and nothing are return
"""

function parse_name_length(newick::String)
    newick = strip(newick)
    print(newick)
    if length(newick) < 1
        println("i an here")
        println(length(newick))
        println(newick)
        return "no_name", 0.0
    end # if length
    if occursin(':',newick)
        name, len = split(newick,':')
        return string(name), parse(Float64, len)
    end # if occusrsin
    return newick, 0.0
end # function

# function make_node(newick::String)
#     node_list = []
#     parts = split(newick, ')')
#     if length(parts) ==1
#         children, current_level = [],parts[1]
#         current_level.split(',')
#         for x in current_level
#             name, length = parse_name_length(x)
#             x = Node()
#             x.name =
#     else
#         println("so why are we here again?")
#
    # if len(parts)
    # len_minus_one = size(parts,1)-1
    # children = list(parse_siblings(join(parts[len_minus_one],')')[2:size(parts,1)]))
    # label = parts[len_minus_one]
    # end #if
    # name, inc_length = parse_name_length(label)
    # parent = Node(name,ones(3,3),missing,children,length(children),false,inc_length,"0",1,0.5,nothing,nothing,true)
    # for x in children
    #     x.mother=parent
    # end #for
    # return parent
#end #function



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
