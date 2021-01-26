# output


```@index
Pages = ["output.md"]
```


## Chains
```@docs
Main.MCPhylo.Chains
```
Construct a `Chains` object that stores MCMC sampler output. Returns an object of type `Chains`.

#### Arguments 
* `iters`: total number of iterations in each sampler run, of which `length(start:thin:iters)` outputted iterations will be stored in the object.

* `params`: number of parameters to store.

* `value`: array whose first, second (optional), and third (optional) dimensions index outputted iterations, parameter elements, and runs of an MCMC sampler, respectively.

* `start`: number of the first iteration to be stored.

* `thin`: number of steps between consecutive iterations to be stored.

* `chains`: number of simulation runs for which to store output, or indices to the runs (default: 1, 2, …).

* `names`: names to assign to the parameter elements (default: `"Param1"`, `"Param2"`, …).




