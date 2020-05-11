using Serialization


"""
There are two newick trees and two julia objects containing these respective trees in
the current data format.

Newick:
tree1.nwk
tree2.nwk

Julia objects
tree1.jls
tree2.jls

You can read in the julia objects as follows:
"""
tree = deserialize("tree1.jls")

 # please supply the correct name of the function
parsed_tree = newick_parser("tree1.nwk")

"""
You can use the deserialized Julia objects to compare them against your trees.

There is a function called 'RF' which you can call on two trees. If the function
returns 0 the trees have an identical topology.
"""
topological_distance = RF(tree, parsed_tree)

"""
If you properly numbered all of the nodes of the tree you parsed from the newick
file, the function 'get_branchlength_vector' returns a vector of all branch lengths
in the tree. You can compare your branch lengths to the original branch lenghts.
Be aware, that the numbering of the nodes might be different and thus the order
of the two vectors. If you abstract from the ordering difference it mights still be
a reasonable indicator.
"""
original_branch = get_branchlength_vector(tree)
parsed_branch = get_branchlength_vector(parsed_tree)
