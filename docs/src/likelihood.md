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
MCPhylo.TreeDistribution
```

### Branch length distributions

```@docs
MCPhylo.CompoundDirichlet
MCPhylo.exponentialBL
MCPhylo.UniformBranchLength
MCPhylo.BirthDeath
MCPhylo.BirthDeathSimplified
```

### Topology distributions

```@docs
MCPhylo.UniformConstrained
MCPhylo.UniformTopology
```

## Gamma Rates

```@docs
MCPhylo.discrete_gamma_rates
```

## Substitution Models

```@docs
MCPhylo.JC
MCPhylo.Restriction
MCPhylo.freeK
MCPhylo.GTR
```
