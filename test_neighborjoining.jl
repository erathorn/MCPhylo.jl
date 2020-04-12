using Serialization


"""
There are three distance matrices and lists of leaves, which you can read into
a julia object as shown here
"""
distance_matrix = deserialize("dm_1.jls")
leaves_list = deserialize("leaves_1.jls")


"""
After creating a tree with the neighbor joining method, you can save the tree as
a newick string and save it into a file as follows
"""
tree = neighbor_joining(distance_matrix, leaves_list)
f = open("newick_output.nwk", "w")
println(f, newick(tree))
close(f)


"""
Using the dendropy python package (https://dendropy.org/) you can compare your tree
to the gold standard one in newick1.nwk (respectively for newick2, and newick3)

Apparently they can also compute a neighbor joining tree from a distance matrix,
if you find it necessary, you can also use these methods for comparison.
"""
