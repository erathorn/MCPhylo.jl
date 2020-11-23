# Tree Functionality
Functions found in all files of Tree folder, ordered according to file of origin.

```@index
Pages = ["Tree.md"]
```

## Converter.jl
```@docs
Main.MCPhylo.from_df
Main.MCPhylo.newick(::T) where T<:Main.MCPhylo.GeneralNode
Main.MCPhylo.to_covariance
Main.MCPhylo.to_covariance_ultra
Main.MCPhylo.to_df
Main.MCPhylo.to_distance_matrix
```

## Node_Type.jl
```@autodocs
Modules = [MCPhylo]
Pages   = ["Node_Type.jl"]
Filter =
```
## Tree_Basics.jl
```@docs
Main.MCPhylo.add_child!
Main.MCPhylo.check_binary
Main.MCPhylo.create_tree_from_leaves
Main.MCPhylo.delete_node!
Main.MCPhylo.force_ultrametric!
Main.MCPhylo.get_branchlength_vector
Main.MCPhylo.get_mother
Main.MCPhylo.get_sister
Main.MCPhylo.get_sum_seperate_length!
Main.MCPhylo.insert_node!
Main.MCPhylo.node_height
Main.MCPhylo.number_nodes!
Main.MCPhylo.path_length
Main.MCPhylo.random_node
Main.MCPhylo.remove_child!
Main.MCPhylo.set_binary!
Main.MCPhylo.set_branchlength_vector!
Main.MCPhylo.tree_height
Main.MCPhylo.tree_length(::T) where T<:Main.MCPhylo.GeneralNode
```

## Tree_Clustering.jl
```@docs
Main.MCPhylo.neighbor_joining
Main.MCPhylo.upgma
```

```
```@autodocs
Modules = [MCPhylo]
Pages   = ["Tree_Clustering.jl"]
Filter =
```
```
## Tree_Consensus.jl
```@autodocs
Modules = [MCPhylo]
Pages   = ["Tree_Consensus.jl"]
Filter =
```
## Tree_Distance.jl
```@autodocs
Modules = [MCPhylo]
Pages   = ["Tree_Distance.jl"]
Filter =
```
## Tree_Ladderizing.jl
```@autodocs
Modules = [MCPhylo]
Pages   = ["Tree_Ladderizing.jl"]
Filter =
```
## Tree_Legacy.jl
```@autodocs
Modules = [MCPhylo]
Pages   = ["Tree_Legacy.jl"]
Filter =
```
## Tree_moves.jl
```@docs
Main.MCPhylo.NNI!
Main.MCPhylo.NNI
Main.MCPhylo.change_edge_length!
Main.MCPhylo.move!
Main.MCPhylo.perform_spr
Main.MCPhylo.randomize!
Main.MCPhylo.risky_SPR!
Main.MCPhylo.risky_SPR
Main.MCPhylo.slide!
Main.MCPhylo.slide
Main.MCPhylo.swing!
Main.MCPhylo.swing
```

## Tree_Pruning.jl
```@docs
Main.MCPhylo.prune_tree!
Main.MCPhylo.prune_tree
```

## Tree_Search.jl
```@docs
Main.MCPhylo.find_binary
Main.MCPhylo.find_num
Main.MCPhylo.find_root
```

## Tree_Traversal.jl
```@docs
Main.MCPhylo.get_leaves
Main.MCPhylo.level_order
Main.MCPhylo.level_traverse
Main.MCPhylo.post_order
Main.MCPhylo.pre_order
```
