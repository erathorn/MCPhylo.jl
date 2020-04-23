
#include("../Tree/Tree_Traversal.jl")
include("../MCPhylo.jl")
using ..MCPhylo

#include("./MCPhylo/src/MCPhylo.jl")
#using .MCPhylo #  ==> depending on the level of embeding

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
    node = newick_parse(content,nothing)


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

# the grand line between functions which are we more certain in and experimental functions on which we are still working on

###########

# TODO: meaningful comment
"""
    newick_parse(newick::String, current_node::Any)

This function parses the string in newick format and returns Node
"""

#parses the string, creating nodes/recursing as needed
function newick_parse(newick::String, current_node::Any,my_nodes::Array{Any,1}=[]) #TODO: INSTEAD OF DOING RECURSION ON ALL DESCENDANTS AT ONCE, RECURSE FOR EVERY CHILD
    #first time function runs on a given string, no node(or an empty node) should be input
    #if no node is input, this happens

    if current_node == nothing
        newick = strip(newick)
        current_node = Node()
        current_node.name = "Root"
        my_nodes::Array{Any,1} = []
    end #if
    println("it begins")
    println("current newick: ",newick)

    while true
        #happens when input string is fully parsed; ends recursion as well
         if newick == "" || newick == ";"
             println("grats, it's over and it went perfectly without any need for further work")
             return current_node, my_nodes
         end #if
         #commas are basically ignored; if there's a parenthesis it means we need to recurse, if there's
         #letters it means we just make a leaf node. This just removes commas from the string.
         if newick[1] == ',' #we don't need these i don't think, just gotta look at the next thing
             newick = SubString(newick,2)
             println("what follows is a sibling")
             println("sibling: ",newick)
         end #if
         if newick[1] == '('
             #this is recursion; if we're looking at an internal node this should happen
             println("ALERT: recursion needed")
             childs_section = match(r"\(([^()]|(?R))*\)",newick) #should return only the descendants of the current child node, check https://regex101.com/r/lF0fI1/1 for proof of this regex working the way it should
             childs_section = childs_section.match
             index=findlast(")",childs_section)[1]+1
             the_rest = SubString(childs_section,index)
             newick = the_rest

             """
             parse_newick((A,B),C)
             (A,B)C;
             add_child!(C, parse_newick(A))
             add_child!(C, parse_newick(B))

             function parse_newick(me)
                 if I am leave
                     me, []
                     just return me
                 else
                     me, my_children = parse_newick(me)
                     if my_children is empty
                         return me, []
                    else
                        add my_children to me
                        return me, []
        test on these:
            A; => Base case A
            (A,B)C;
            ((A,B),C)E;
            (C,(A,B))E;
             ((A,B),(C,D))E;

        test on these ^
        newick(rootnode)

        (A,B)C
        (B,A)C
        [n.name for n in post_order!]

        UI Level
            1 sanity check => is there a semicolon
            2 root = internal recursion
            ####
            3 set_binary!(root)
            4 number_nodes!(root)

            return root


             add_child!(E, (A,B))
             add_child!(E, (C,D))

             ((A,B),(C,D));

             using Juno
             @enter newick_parser()


             """
             childs_section = SubString(childs_section,2,findlast(')',childs_section)-1)
             push!(my_nodes,current_node)
             cur_child = Node()
             add_child!(current_node,newick_parse(string(childs_section),cur_child,my_nodes)[1]) #where the magic happens
             if newick == ""
                 return current_node,my_nodes
             end #if
                if newick[1] == ':' #parses name and length of current node
                    node_boarder = match(r"[();,]",newick).offset
                    name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
                    current_node.inc_length = length
                     push!(my_nodes,current_node)
                    newick = SubString(newick,node_boarder)
                end #if
            end #if

            if occursin(r"^[0-9A-Za-z_|]+",string(newick[1])) || newick[1] == ':'
                println("we found a leaf")

                #this should happen if the current node's a leaf node
                if match(r"[();,]",newick) == nothing #if we're parsing the last leaf node; prevents index out of bounds errors
                    println("leaf is: ", newick)
                    endpoint = lastindex(newick)
                    name,length = parse_name_length(string(SubString(newick,1,endpoint)))
                    cur_child = Node()
                    cur_child.name = name
                    cur_child.inc_length = length
                    add_child!(current_node,cur_child)
                    push!(my_nodes,cur_child)
                    newick = ""
                else
                    node_boarder = match(r"[();,]",newick).offset #handles creation of every leaf node BUT the last one, which is handled ^
                    println("leaf is: ", string(SubString(newick,1,node_boarder-1)))
                    name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
                    cur_child = Node()
                    cur_child.name = name
                    cur_child.inc_length = length
                    add_child!(current_node,cur_child)
                    push!(my_nodes,cur_child)
                    newick = SubString(newick,node_boarder)
            end #if
        end #while
    end #function
end


"""
    parse_name_length(newick::String)

This function parses two optional elements of the tree, name and length. In case, when neither of this is provided, empty string and nothing are return
"""

function parse_name_length(newick::String)
    newick = strip(newick)
    #name, my_length = split(newick, ":")
    # if name == "" #=> no name is given, use dummy

    # if my_length == "" #=> no length is given, use dummy
    if length(newick) < 1
        return "nameless", 1.0
    end # if length
    if occursin(':',newick)
        name, len = split(newick,':')
        return string(name), parse(Float64, len)
    end # if occusrsin
    return newick, 1.0
end # function

###
# TESTING THE FUNCTIONS


the_string = "(Swedish_0:0.1034804,(Welsh_N_0:0.1422432,(Sardinian_N_0:0.02234697,(Italian_0:0.01580386,Rumanian_List_0:0.03388825):0.008238525):0.07314805):0.03669193,(((Marathi_0:0.04934081,Oriya_0:0.02689862):0.1193376,Pashto_0:0.1930713):0.05037896,Slovenian_0:0.0789572):0.03256979);"
the_string_2 = "((((Dra.NORTHERN_DRAVIDIAN.KURUKH:0.823074875737579,Dra.SOUTH_CENTRAL_DRAVIDIAN.DORLI_GONDI:0.28822305943964277)51:0.9096574009527566,Dra.SOUTH_CENTRAL_DRAVIDIAN.KUVI:0.15108563362392177)64:0.6885640496300766,(Dra.SOUTH_CENTRAL_DRAVIDIAN.SOUTH_EASTERN_GONDI:0.9621464574054764,((Dra.SOUTH_CENTRAL_DRAVIDIAN.KUI:0.7660539461805399,(Dra.SOUTHERN_DRAVIDIAN.TAMIL:0.22241296989623094,(Dra.SOUTH_CENTRAL_DRAVIDIAN.WESTERN_GONDI:0.700208475875382,Dra.SOUTH_CENTRAL_DRAVIDIAN.ADILABAD_GONDI:0.00636902086929546)45:0.2568376168446603)46:0.31195689986105074)50:0.5280782014653551,(Dra.SOUTH_CENTRAL_DRAVIDIAN.PENGO:0.3805497140792932,Dra.SOUTHERN_DRAVIDIAN.IRULA:0.24313087586739268)47:0.5514806675227825)54:0.40644694508297363)56:0.7481616281899098)70:0.6308545527775923,(((((Dra.SOUTH_CENTRAL_DRAVIDIAN.SOUTHERN_GONDI:0.14981318626586068,Dra.SOUTHERN_DRAVIDIAN.TULU:0.3189677711581785)52:0.18701136213491754,Dra.SOUTH_CENTRAL_DRAVIDIAN.ABUJHMARIA:0.9570448784020357)55:0.7996910080326368,(Dra.CENTRAL_DRAVIDIAN.NORTHWESTERN_KOLAMI:0.46552096552619276,Dra.SOUTHERN_DRAVIDIAN.PALIYAN:0.005423542225934379)44:0.5563938421755457)68:0.6220646419257747,((Dra.NORTHERN_DRAVIDIAN.KURUX_NEPALI:0.4546646140265832,Dra.SOUTHERN_DRAVIDIAN.RAVULA:0.9779586706039829)60:0.7619109389692604,((Dra.NORTHERN_DRAVIDIAN.KUMARBHAG_PAHARIA:0.8299924660663025,Dra.SOUTH_CENTRAL_DRAVIDIAN.TELUGU_2:0.6170375335378135)49:0.4061078279052086,(Dra.SOUTHERN_DRAVIDIAN.OLD_TAMIL:0.25825955834703956,Dra.SOUTH_CENTRAL_DRAVIDIAN.NORTH_BASTAR_GONDI:0.6154523803114922)59:0.3592581704448605)62:0.20575614618578816)67:0.8445878101331729)69:0.4277614887197028,((Dra.SOUTH_CENTRAL_DRAVIDIAN.TELUGU:0.05936684953610001,(Dra.SOUTHERN_DRAVIDIAN.KURUMBA_ALU:0.39761369422055753,Dra.SOUTH_CENTRAL_DRAVIDIAN.NORTHERN_GONDI:0.3512494180801804)43:0.024384381704830566)61:0.6159956440241587,((((Dra.SOUTHERN_DRAVIDIAN.BETTA_KURUMBA:0.556468834623374,((Dra.NORTHERN_DRAVIDIAN.SAURIA_PAHARIA:0.7344814114120296,Dra.SOUTH_CENTRAL_DRAVIDIAN.SOUTHERN_GONDI_2:0.7759235873330327)53:0.14971493972622438,Dra.SOUTHERN_DRAVIDIAN.KORAGA_KORRA:0.7802715667203111)57:0.5978920373858239)58:0.5125071989391269,(Dra.SOUTH_CENTRAL_DRAVIDIAN.SOUTH_BASTAR_GONDI:0.4634536648108183,((Dra.NORTHERN_DRAVIDIAN.BRAHUI:0.4481324381951019,Dra.SOUTHERN_DRAVIDIAN.KANNADA:0.4493635834393417)40:0.5017018479452859,Dra.SOUTHERN_DRAVIDIAN.MALAYALAM:0.6569175054247594)48:0.36416624962676897)63:0.9195968543725148)65:0.31970982925608965,(Dra.SOUTHERN_DRAVIDIAN.TODA:0.8335870971561289,(Dra.CENTRAL_DRAVIDIAN.PARJI:0.05328619625991444,Dra.SOUTHERN_DRAVIDIAN.ULLATAN:0.9638487161160924)39:0.6271628983114025)41:0.19457266408948884)66:0.4234684186473759,((Dra.SOUTH_CENTRAL_DRAVIDIAN.MANDA_INDIA:0.8101005415604958,Dra.SOUTHERN_DRAVIDIAN.BADAGA:0.13465908110800426)42:0.8025801334975355,Dra.SOUTH_CENTRAL_DRAVIDIAN.HILL_MARIA_GONDI:0.8975962746953261)71:0.26228038217619115)72:0.23390480938782518)73:0.22512880787981576)74:0.7485699845253884)75:0.5;"
println("begin of new test is here")
println("")
println("")
println("")
println("")
println("")
println("")

# -> dendroscpoe

the_easiest_one = "(A,B)C;"
# (,,(,));                               no nodes are named
# (A,B,(C,D));                           leaf nodes are named
named_leaf_nodes = "(A,B,(C,D));"           #makes nodelist,makes correct tree, names leaf nodes correctly
# (A,B,(C,D)E)F;                         all nodes are named
all_nodes_named = "(A,B,(C,D)E)F;"         #makes correct tree, names leaf nodes correctly, makes nodelist correctly, does NOT name internal nodes correctly
# (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
almost_all_distances = "(:0.1,:0.2,(:0.3,:0.4):0.5);" #makes correct tree, makes nodelist correctly
# (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
all_distance_to_parent = "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"#nodelist, tree, leaf names all fine
# (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
distance_and_leaf_name = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);" #nodelist, tree, leaf names all fine
# (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
distance_and_all_names = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"#nodelist, tree, leaf names all fine
# ((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;    a tree rooted on a leaf node (rare)
the_rare_one = "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;" #nodelist, tree, leaf names all fine

#newick(rootnode) <--- will return a newick string representation of a tree; can use to test instead of trying to print every node
#post_order <--- will return a vector of nodes in post order
#pre_order <--- same, but pre order
#level_order <--- same, but level order

test1,nodelist = newick_parse(the_string_2,nothing)
println("the useless root node is: ", test1)
println("LIST FOLLOWS")
println("")
println("")
println("")
println("")
println("")
println("")
println("")
println(nodelist)
 kids = test1.children
  for x in kids
     println("child of the above node is: ",x)
     grandkids = x.children
     for y in grandkids
         println("children of ",x.name, " are :", y)
         greatgrandkids = y.children
         for z in greatgrandkids
             println("children of ",y.name, " are: ",z)
             greatgreat = z.children
             for blah in greatgreat
                 println("children of ",z.name, " are: ",blah)
                 greatgreatgreat = blah.children
                 for blahblah in greatgreatgreat
                     println("children of ",blah.name, " are: ",blahblah)
                 end
             end
         end
     end
 end


###
# SEMMINGLY USELESS FUNCTIONS FOR NOW
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
#
#
# function parsing_the_newick(newick::String,current_node::Any,count::Integer)
#
#     if current_node == nothing
#         count = 0
#         current_node = Node()
#         println()
#         println()
#         println()
#     end # setting things up
#
#     # the child case
#
#     if newick[1] == '(' || newick[1] == ','
#
#         newick = SubString(newick,2)
#
#         if newick[1] == '('
#             # YOUR RECURSION IS HERE
#             child = parsing_the_newick(string(newick),current_node,count)
#             add_child!(current_node,child)
#         end # if (the recursive call one)
#
#         node_boarder = match(r"[();,]",newick).offset
#         println("Trying to detect the name and length from the ", string(SubString(newick,1,node_boarder-1)))
#         name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
#
#         child = Node()
#         child.name = name
#         child.inc_length = length
#         child.num = count
#         count+=1
#         add_child!(current_node,child)
#         newick = SubString(newick,node_boarder)
# end # if "("
#
# if newick[1] == ')'
#         newick = string(SubString(newick,2))
#         println("Continiue to parse... ", newick)
# end # if with the ")" bracket
#
# if occursin(r"^[0-9A-Za-z_|]+",string(newick[1])) || newick[1] == ':'
#         println("We're getting the information of the current node")
#         node_boarder = match(r"[();,]",newick).offset
#         println("Trying to detect the name and length from the ", string(SubString(newick,1,node_boarder-1)))
#         name,length = parse_name_length(string(SubString(newick,1,node_boarder-1)))
#         current_node.name = name
#         current_node.inc_length = length
#         current_node.num = count
#         count+=1
#         newick = SubString(newick,node_boarder)
#         println("Information about the current node was written down. Continue to parse... ",newick)
# end # if length one
#     if newick[1] == ';'
#         println("We have reached the end!")
#     end # the last one
# return newick, current_node
# end # the function
# function parse_siblings(newick::String)
#     bracket_lvl = 0
#     current = []
#     for c in  (newick * ',')
#         if c == ','
#             if bracket_lvl == 0
#                 yield(make_node(join(current,""))) => this is not Julia
#                 current = []
#             else
#                 if c == '('
#                     bracket_lvl += 1
#                 elseif c == ')'
#                     bracket_lvl -= 1
#                 end #elseif
#             end #if
#             push!(current,c)
#         end #if
#     end #for
# end #function
