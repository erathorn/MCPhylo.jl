# Model

The Model type is designed to store the set of all model nodes, including parameter set \Theta. In particular, it stores Dependent type objects in its nodes dictionary field. Valid models are ones whose nodes form directed acyclic graphs (DAGs). Sampling functions ``\{f_j\}_{j=1}^{B}`` are saved as Sampler objects in the vector of field samplers. Vector elements j=1,\ldots,B correspond to sampling blocks ``\{\Theta_j\}_{j=1}^{B}``.

## Dependent

```@docs
MCPhylo.Logical
MCPhylo.Stochastic
```

## Graph

```@docs
MCPhylo.draw
MCPhylo.graph2dot
```

## Initialization

```@autodocs
Modules = [MCPhylo]
Pages   = ["initialization.jl"]
Filter =
```

## MCMC

```@docs
Main.MCPhylo.mcmc
```

## Model Building

```@docs
MCPhylo.Model
```
