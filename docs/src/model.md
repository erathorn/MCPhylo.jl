# Model

The Model type is designed to store the set of all model nodes, including parameter set \Theta. In particular, it stores Dependent type objects in its nodes dictionary field. Valid models are ones whose nodes form directed acyclic graphs (DAGs). Sampling functions ``\{f_j\}_{j=1}^{B}`` are saved as Sampler objects in the vector of field samplers. Vector elements j=1,\ldots,B correspond to sampling blocks ``\{\Theta_j\}_{j=1}^{B}``.

## Dependent

```@docs
MCPhylo.Logical
MCPhylo.Stochastic
MCPhylo.setinits!
MCPhylo.setmonitor!
MCPhylo.update!
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

## Model

```@docs
Main.MCPhylo.Model
Main.MCPhylo.Base.keys
Main.MCPhylo.Base.show
Main.MCPhylo.showall
```

## Simulation
```@docs
Main.MCPhylo.gettune
Main.MCPhylo.settune!

Main.MCPhylo.gradlogpdf
Main.MCPhylo.gradlogpdf!
Main.MCPhylo.logpdf
Main.MCPhylo.logpdf!
Main.MCPhylo.sample!(::Main.MCPhylo.Model, ::Integer)
Main.MCPhylo.unlist
Main.MCPhylo.relist
Main.MCPhylo.relist!
Main.MCPhylo.update!
```
