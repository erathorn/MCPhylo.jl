# model
```@index
Pages = ["model.md"]
```

## Model

```@docs
Main.MCPhylo.Model
```

#### Arguments

* `iter`: current iteration of the MCMC simulation.

* `burnin`: number of initial draws to be discarded as a burn-in sequence to allow for convergence.

* `samplers`: block-specific sampling functions.

* `nodes...`: arbitrary number of user-specified arguments defining logical and stochastic nodes in the model. Argument values must be `Logical` or `Stochastic` type objects. Their names in the model will be taken from the argument names.


[Back to top.](@ref model)