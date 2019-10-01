# MCPhylo

This package does phylogenetic computations in Julia. It comes with its own tree
module which currently supports the output of newick strings for tree structures.
The goal is to facilitate phylogenetic computations in computational historical
linguistics. To facilitate the inference of phylogenetic trees,
_Probabilistic Path Hamiltonian Dynamics_ (https://arxiv.org/pdf/1702.07814.pdf)
are implemented.

**This package is currently under heavy development.**

**It is written in Julia 1.2.0. This package is not backwards compatible!**

**This module is build on a forked instance of Mamba 0.12.0**

## How to

**General Information**

A proper installation of this package is not supported yet. Please install all the
necessary packages listed in `REQUIRE`. Next, download this
package and place it in your workspace. Then use it as shown in `tester.jl`.

An initial running version is shown in the file `tester.jl`.
The general setup is as in the original Mamba package.

In addition to the distributions supported by `Mamba` and the `Distributions` package
there are currently three distributions for trees supported by the `MCPhylo` package.

* CompoundDirichlet (Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030))
* Strict Molecular Clock -Birth Death following Yang & Rannala 1997 (doi.org/10.1093/oxfordjournals.molbev.a025811)
* Strict Molecular Clock - Simplified Birth Death Yang & Rannala 1996 (doi.org/10.1007/BF02338839)

___Caveat:___ The Molecular Clock Priors are not yet tested. There is a uniform prior on tree topologies.

Only Nexus files are with binarised cognate judgments are currently supported.

### Recipe

_Note:_ Only the steps which are different from the standard `Mamba` model setup
are listed here.

1. Read in the nexus file using the `make_tree_with_data` function. It returns a
tree object including the cognate information.
2. Specify the Model using the standard `Mamba` syntax for specifying a model. To properly
use the Probabilistic Path approach, the objective should be a `PhyloDist` distribution object.
3. The Probabilistic Path sampler is selected using the `ProbPathHMC` sampler. The parameters of this
function are in this order:
    1. The symbol specifying the sampling target
    2. Number of leap-prog steps
    3. leap-prog stepsize
    4. smoothing threshold
    5. the `:provided` symbol to use the gradient which is provided by the package
_Note:_ The Probabilistic Path sampler will not work, if you do not use the provided Gradient function
4. The standard `mcmc` function from `Mamba` takes an extra boolean argument `trees` indicating if the
sampled trees should be stored. If set to _true_ the trees will be stored. The default is _false_.
5. You can flush the model parameters and the sampled trees to a file, using the `to_file` function.
It takes as a first argument an MCMC object and second a path to a folder where the results should be stored.
