module Tree_Module
include("Tree_Basics.jl")
include("Converter.jl")

using .Tree_Basics
using .Converter

export create_tree_from_leaves, post_order, tree_length, tree_height,
       path_length, get_leaves, to_df, Node, add_child!


end # Tree_Module
