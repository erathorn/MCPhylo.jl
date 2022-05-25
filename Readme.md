# MCPhylo

[![CI](https://github.com/erathorn/MCPhylo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/erathorn/MCPhylo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/erathorn/MCPhylo.jl/branch/main/graph/badge.svg?token=05DQRGQAIK)](https://codecov.io/gh/erathorn/MCPhylo.jl)

This package offers the functionality for phylogenetic computations in Julia. It is an extension of the `Mamba` package which implements different variants of Markov Chain Monte Carlo (MCMC) sampling for Bayesian analysis. ([https://mambajl.readthedocs.io/en/latest/](https://mambajl.readthedocs.io/en/latest/))
**MCPhylo** extends `Mamba` to perform phylogenetic computations.

**This module is build on a forked instance of Mamba 0.12.0**

The documentation can be found [here](https://erathorn.github.io/MCPhylo.jl/dev/).

------

## Citation

If you use this software in an academic publication please cite the following paper
No-U-Turn sampling for phylogenetic trees, Wahle (2021).  ([bioRxiv Paper](https://doi.org/10.1101/2021.03.16.435623))

------

## General Information

This package can be installed via the official Julia PackageManager.

```julia
using Pkg
Pkg.add("MCPhylo")
```

The latest functionalities can be obtained by installing the it directly from the GitHub source.

```julia
using Pkg
Pkg.add("https://github.com/erathorn/MCPhylo.jl")
```
