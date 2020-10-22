# Tree Functionalities

There are several tree functionalities already implemented.

## First order tree Functions

- tree_length!
- tree_height!
- node_height!
- node_distance
- get_path
- get_sister!
- get_mother!
- find_lca

## Tree Manipulation

- NNI! - do a random nearest neighbor interchange move on the tree
- slide! - do a slide move on a random node
- swing! - do a swing move on a random node
- change_edge_length! - change the edge length of a random node
- ladderize_tree - return a ladderized version of the tree

## Tree Traversal

- post_order!
- pre_order!
- level_order!

## Tree Distance

- RF - Calculate the Robinson Foulds distance between two trees
- BHV_bounds - calculate lower and upper bound on the geodesic in the Billera-Holmes-Vogtman space

## Tree Construction

- upgma - Construct a tree from a distance matrix using the UPGMA algorithm
- neighbor_joining -  Construct a tree from a distance matrix using the neighbor joining algorithm

## Tree Conversion

- to_df - return a matrix representation of the tree structure
- from_df - return a tree structure from a matrix representation
- to_covariance - Calcualte the variance-covariance matrix from a tree
- to_distance_matrix - Calculate the distance matrix over the set of leaves.
- newick - return a newick representation of the tree

## Data Parsing

### Wrapper

- make_tree_with_data - create a random tree on the data defined in file. Can either be a Nexus or a CSV file

### Non Wrapper

- ParseNexus - Parses a NEXUS file which stores the data. The file should follow the conventions used for MrBayes.
- ParseCSV - Parses a CSV file which stores the data.
- ParseNewick - Parse a file, containing Newick strings.
- parsing_newick_string -  parse a newick string into a tree (*Disclaimer: This function does not check the validity of the newick string*)

## Distributions on Trees

- CompoundDirichlet (Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030))