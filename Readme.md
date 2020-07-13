# MCPhylo

This package does phylogenetic computations in Julia. It comes with its own tree
module which currently supports the output of newick strings for tree structures.
The goal is to facilitate phylogenetic computations in computational historical
linguistics. To facilitate the inference of phylogenetic trees,
_Probabilistic Path Hamiltonian Dynamics_ (https://arxiv.org/pdf/1702.07814.pdf)
are implemented.

**This package is currently under heavy development.**

**It is written in Julia 1.3.1. This package is not backwards compatible!**

**This module is build on a forked instance of Mamba 0.12.0**

## How to

**General Information**

Installation of the package should work.

```julia
using Pkg
Pkg.add("https://githubg.com/erathorn/JuliaTree")
```


**Note** This package uses multithreading. (https://docs.julialang.org/en/v1/base/multi-threading/) and is thus not backwards compatible.

An initial running version is shown in the file `PNUTSExample.jl`. The general setup is as in the original Mamba package.

In addition to the distributions supported by `Mamba` and the `Distributions` package there are currently three distributions for trees supported by the `MCPhylo` package.

* CompoundDirichlet (Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030))
* Strict Molecular Clock -Birth Death following Yang & Rannala 1997 (doi.org/10.1093/oxfordjournals.molbev.a025811)
* Strict Molecular Clock - Simplified Birth Death Yang & Rannala 1996 (doi.org/10.1007/BF02338839)

___Caveat:___ The Molecular Clock Priors are not yet tested. There is a uniform prior on tree topologies.

Nexus and CSV files with binarized cognate data are supported.

### Recipe

The standard `mcmc` function from `Mamba` takes an extra boolean argument `trees` indicating if the sampled trees should be stored. If set to _true_ the trees will be stored. The default is _false_.

You can flush the model parameters and the sampled trees to a file, using the `to_file` function. It takes as a first argument an MCMC object and second a path to a folder where the results should be stored. This file can be read by the Tracer softward (https://github.com/beast-dev/tracer/)
