# Likelihood
There are several prior distributions for tree structures implemented in this package.
Also there is a function which calculates the likelihood of a tree given a model using
Felsensteins algorithm.

## Likelihood Calculator functionality
```@autodocs
Modules = [MCPhylo]
Pages   = ["Likelihood/LikelihoodCalculator_Node.jl"]
Filter =
```

## Prior
```@docs
Main.MCPhylo.CompoundDirichlet
Main.MCPhylo.exponentialBL
```
