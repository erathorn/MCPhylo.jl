module Tree_Module
include("Tree_Basics.jl")
include("Converter.jl")
include("Tree_moves.jl")

using .Tree_Basics
using .Converter
using .Tree_moves

export create_tree_from_leaves, post_order, tree_length, tree_height,
       path_length, get_leaves, to_df, Node, add_child!, remove_child!, NNI!,
       slide!, swing!, to_newick


end # Tree_Module
