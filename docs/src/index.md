```@contents
```

# Adding MCPhylo
```@repl
using Pkg
Pkg.add("MCPhylo")
```

# Parsing Functionality
```@docs
MCPhylo.ParseNewick
```
```@docs
MCPhylo.ParseCSV
```
```@docs
MCPhylo.ParseNexus
```

# Tree Manipulation
```@docs
MCPhylo.move!
```
```@docs
MCPhylo.merge_trees!
```
```@docs
MCPhylo.ladderize_tree!
```
```@docs
MCPhylo.prune_tree!
```
```@docs
MCPhylo.prune_tree
```
```@docs
MCPhylo.add_child!
```
```@docs
MCPhylo.remove_child!
```
```@docs
MCPhylo.create_tree_from_leaves
```
```@docs
MCPhylo.create_tree_from_leaves_bin
```
```@docs
MCPhylo.insert_node!
```
```@docs
MCPhylo.number_nodes!
```
```@docs
MCPhylo.perform_spr
```
```@docs
MCPhylo.risky_SPR
```
```@docs
MCPhylo.risky_SPR!
```

# Return Functions
```@docs
MCPhylo.newick
```
```@docs
MCPhylo.tree_height
```
```@docs
MCPhylo.level_order
```
```@docs
MCPhylo.tree_length
```
```@docs
MCPhylo.get_leaves
```
```@docs
MCPhylo.get_bipartitions
```
```@docs
MCPhylo.get_sister
```
```@docs
MCPhylo.get_branchlength_vector
```
```@docs
MCPhylo.get_mother
```
```@docs
MCPhylo.get_sum_seperate_length!
```