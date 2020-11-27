# MCPhylo

This package does phylogenetic computations in Julia. It is an extension of the `Mamba` package which does Markov Chain Monte Carlo (MCMC) sampling for Bayesian analysis. ([https://mambajl.readthedocs.io/en/latest/](https://mambajl.readthedocs.io/en/latest/))
**MCPhylo** extends `Mamba` by a tree module to perform phylogenetic computations.

which currently supports the output of newick strings for tree structures.
The goal is to facilitate phylogenetic computations in computational historical
linguistics. To facilitate the inference of phylogenetic trees,
_Probabilistic Path Hamiltonian Dynamics_ (https://arxiv.org/pdf/1702.07814.pdf)
are implemented.

**This package is currently under heavy development.**

**This package needs at least Julia 1.3.1. This package is not backwards compatible!**

**This module is build on a forked instance of Mamba 0.12.0**

**Note** This package uses multithreading. (https://docs.julialang.org/en/v1/base/multi-threading/)



## General Information

In order to use the current version of the package clone the repo and place it into your current working directory.

```julia
include("./src/MCPhylo.jl")
using .MCPhylo
```

Installation of the package may also work.

```julia
using Pkg
Pkg.add("https://github.com/erathorn/JuliaTree")
```


The setup of a model is as in the original Mamba package.

Nexus and CSV files with binarized cognate data are supported.

## New Functions

The standard `mcmc` function from `Mamba` takes an extra Boolean argument `trees` indicating if the sampled trees should be stored. If set to _true_ the trees will be stored. The default is _false_.

You can flush the model parameters and the sampled trees to a file, using the `to_file` function. It takes as a first argument an MCMC object and second a path to a folder where the results should be stored. This file can be read by the Tracer software (https://github.com/beast-dev/tracer/). Additionally, if trees are stored it will create a file with newick strings of these trees.

### New & Adjusted Samplers

`PNUTS ` is a sampler which does *Phylogenetic No-U-Turn sampling* (Wahle (forthcomming)). It samples tree stochastic nodes.

`RWM` Random walk metroplis hastings sampling can work with trees now. For numerical nodes the sampler and the function signature is as in the original Mamba package. For tree structures the signature is slightly different: `RWM(:tree, :all)` or `RWM(:tree, [:NNI, :Swing])` The first variant uses all available tree manipulation moves (see [Tree Manipulation](TreeStuff.md#Tree Manipulation)), the second variant only makes use of a user defined subset of these moves. Ladderization of the tree is not an eligible tree manipulation move

`NUTS` can take the argument `dtype=:Zygote` to use Zygote for the calculation of the gradient. The default is finite differencing.

`Slice` can also sample trees. It does a slice sampling operation on the branch lengths of the tree.
## Tree Functionalities

For the available tree functionalities see: [Tree Functionalities](TreeStuff.md), and for exact documentation of the functions involved, see the [Tree Functionality page] (Tree.md) linked along the left-hand side of the page.

```@contents
Pages = ["Tree.md", "Likelihood.md", "Parser.md", "distributions.md", "model.md", "output.md", "Sampler.md", "samplers.md", "Substitution.md", "Utils.md"]
```
