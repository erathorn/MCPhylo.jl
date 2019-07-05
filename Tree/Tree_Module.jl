module Tree_Module


using Markdown
using Random
using Distributions
using DataFrames



include("Tree_Basics.jl")
include("Converter.jl")
include("Tree_moves.jl")

export create_tree_from_leaves, post_order, tree_length, tree_height,
       path_length, get_leaves, to_df, Node, add_child!, remove_child!, NNI!,
       slide!, swing!, to_newick, find_by_name, find_by_binary, find_by_root, get_branchlength_vector!
       get_sum_seperate_length!, set_branchlength_vector!


end # Tree_Module
