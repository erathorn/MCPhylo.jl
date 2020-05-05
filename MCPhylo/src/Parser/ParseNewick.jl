include("../MCPhylo.jl")
using ..MCPhylo

# using .MCPhylo #  ==> depending on the level of embeding



"""
    load_newick(filename::String)

This function loads a newick from file
"""
# rn it assumed that there are no extra line breaks (\n\n) and there is only one tree pro file. That should be fixed.

function load_newick(filename::String)
    """
    let this function always return a list of strings.

    if there is more than one tree in the file, the newick parser can just return
    a list of trees.

    """
    open(filename, "r") do file
        global content = readlines(file)
    end
    content[1]
end

"""
    is_valid_newick_string(newick::String)

This function checks if the given string is valid: is the brackets number matches and if the string ends with ";"
"""

function is_valid_newick_string(newick::String) #TODO: possibly not necessary; could be done as part of recursion, MAYBE
    # Step one: does the stripped string ends with ';'
    if endswith(strip(newick),";")
        # Step two: check for the equal amount of brackets
        bracket_level = 0
        for letter in newick # => char is a type in Julia
            if letter == '('
                bracket_level += 1
            elseif letter == ')'
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


"""
    parse_name_length(newick::String)

This function parses two optional elements of the tree, name and length. In case, when neither of this is provided, empty string and nothing are return
"""

function parse_name_length(newick::String)
    """
    always do split(newick, ":"), there is no error if ":" does not occur, it
    will just return an array containing the newick string

    r = split("", ":")
    r
    > [""]

    This streamlines the function
    """
    if occursin(':',newick)
        name, len = split(strip(newick),':')
        if name == ""
            name = "nameless"
        end # if
        if len == ""
            len = 1.0
        end # if
        if occursin(';',len)
            len = chop(len)
        end # if
    return string(name), parse(Float64, len)
    end # main if

    if length(newick)<1
        return "nameless",1.0
    else
        return string(newick),1.0
    end #if-else
end # function



function parsing_newick_string(newick::String)
    """
    The comming PR of the Neighbor Joining branch contains cleaner node constructors
    use them later on
    """
    if newick[end] == ';' #no need for semicolon
        newick = chop(newick)
    end #if

    if newick[1] != ')' && occursin(r"^[a-zA-Z]+[:]?[0-9]*",newick)
        leaf_node = Node()
        name,len = parse_name_length(newick)
        leaf_node.name = name
        leaf_node.inc_length = len
        return leaf_node
        # base case; only triggered at end of recursion OR if a single node-tree is input

    else
        current_node = Node()
        childrenstring_with_parenthesis = (match(r"\(([^()]|(?R))*\)",newick)).match #returns section of newick corresponding to descendants of current node, check https://regex101.com/r/lF0fI1/1
        index=findlast(')',childrenstring_with_parenthesis)[1]
        childrenstring = SubString(childrenstring_with_parenthesis,2,index-1) #... so that we can remove the superfluous parentheses here


        """
        move this for loop into a separate function
         -> easier maintainance
        """
        child_list = []
        counter = ""
        bracket_depth = 0
        for x in (childrenstring * ",") # splits string identified above into a list, where each element corresponds to a child of current_node
            if x == ',' && bracket_depth == 0
                push!(child_list,counter)
                counter = ""
                continue
            end #if
            if x == '('
                bracket_depth += 1
            end #if
            if x == ')'
                bracket_depth -= 1
            end #if
            counter = counter * x
        end #for

        for x in child_list #recursion happens here
            add_child!(current_node,parsing_newick_string(x))
        end #for

        child_list = []
        info_of_current_node = split(newick,")") #info of current node should always follow the last ")"
        if lastindex(info_of_current_node) == 1
            name,length = parse_name_length(newick)
        else
            test = string(info_of_current_node[lastindex(info_of_current_node)])
            name,len = parse_name_length(string(test))
            current_node.name = name
            current_node.inc_length = len
        end #else

        return current_node
    end #recursion part
    return nothing #if this happens something went wrong, could do exception handling and whatnot eventually, possible #TODO
    """
    exceptions in julia are raised using the "throw()" command.
    """
end #function



"""
    ParseNewick(filename::String)

This function parses a Newick file
"""
function ParseNewick(filename::String)
    content = load_newick(filename)
    println(typeof(content))
    if !is_valid_newick_string(content)
        throw("$filename is not a Newick file!")
    end # if
    newick(parsing_newick_string(string(content)))
end
