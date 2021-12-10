
# MCPhylo

This package does phylogenetic computations in Julia. It is an extension of the `Mamba` package which does Markov Chain Monte Carlo (MCMC) sampling for Bayesian analysis. ([https://mambajl.readthedocs.io/en/latest/](https://mambajl.readthedocs.io/en/latest/))
**MCPhylo** extends `Mamba` by a tree module to perform phylogenetic computations.

which currently supports the output of newick strings for tree structures.
The goal is to facilitate phylogenetic computations in computational historical
linguistics. To facilitate the inference of phylogenetic trees,
_Probabilistic Path Hamiltonian Dynamics_ (<https://arxiv.org/pdf/1702.07814.pdf>)
are implemented.

**This module is build on a forked instance of Mamba 0.12.0**

# General Information

The stable version can  be installed via the Julia package manager.

```julia
using Pkg
Pkg.add("MCPhylo")
```

In order to use the latest version you can install it directly from its GitHub source.

```julia
using Pkg
Pkg.add("https://github.com/erathorn/MCPhylo.jl")
```

## Tree Functionalities

The tree functionalities within this package are reexported from `MCPhyloTree`[https://erathorn.github.io/MCPhyloTree.jl/stable/]. The the respective documentation for more details.
