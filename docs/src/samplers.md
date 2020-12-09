# samplers
Functions found in all files of samplers folder, ordered according to file of origin.


## abc.jl
```@docs
Main.MCPhylo.ABC
```
#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `scale` : scaling value or vector of the same length as the combined elements of nodes `params` for the `proposal` distribution.  Values are relative to the unconstrained parameter space, where candidate draws are generated.

* `summary` : function that takes a vector of observed or simulated data and returns a summary statistic or vector of statistics.

* `epsilon` : target tolerance for determining how similar observed and simulated data summary statistics need to be in order to accept a candidate draw.

* `kernel` : weighting kernel density of type `Biweight`, `Cosine`, `Epanechnikov`, `Normal`, `SymTriangularDist`, `SymUniform`, or `Triweight` to use in measuring similarity between observed and simulated data summary statistics.  Specified `epsilon` determines the standard deviation of Normal kernels and widths of the others.

* `dist` : positive function for the kernel density to compute distance between vectors of observed (`Tobs`) and simulated (`Tsim`) data summary statistics (default: Euclidean distance).

* `proposal` : symmetric distribution of type `Biweight`, `Cosine`, `Epanechnikov`, `Normal`, `SymTriangularDist`, `SymUniform`, or `Triweight` to be centered around current parameter values and used to generate proposal draws.  Specified `scale` determines the standard deviations of Normal proposals and widths of the others.

* `maxdraw` : maximum number of unaccepted candidates to draw in each call of the sampler.  Draws are generated until one is accepted or the maximum is reached.  Larger values increase acceptance rates at the expense of longer runtimes.

* `nsim` : number of data sets to simulate in deciding whether to accept a candidate draw.  Larger values lead to closer approximations of the target distribution at the expense of longer runtimes.

* `decay` : if `0 < decay <= 1`, the rate at which internal tolerances are monotonically decreased from the initial distance between observed and simulated summary statistics toward the maximum of each subsequent distance and `epsilon`; if `decay = 0`, internal tolerances are fixed at `epsilon`.

* `randeps` : whether to perturb internal tolerances by random exponential variates.

* `args...` : additional keyword arguments to be passed to the `dist` function.

## amm.jl
```@docs
Main.MCPhylo.AMM
```

#### Arguments

* `params` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `Sigma` : covariance matrix for the non-adaptive multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.

* `adapt` : type of adaptation phase.  Options are
    * `:all` : adapt proposal during all iterations.
    * `:burnin` : adapt proposal during burn-in iterations.
    * `:none` : no adaptation (multivariate Metropolis sampling with fixed proposal).
	
* `args...` : additional keyword arguments to be passed to the `AMMVariate` constructor.


```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.AMMVariate, ::Function; ::Bool=true)
```

#### Arguments

* `v` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the `v` argument in a successive call to the function will contain the `tune` field returned by the previous call.
	
* `adapt` : whether to adaptively update the proposal distribution.

## amwg.jl
```@docs
Main.MCPhylo.AMWG
```

#### Arguments

* `params`:  stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `sigma`: scaling value or vector of the same length as the combined elements of nodes 

* `params`, defining initial standard deviations for univariate normal proposal distributions. Standard deviations are relative to the unconstrained parameter space, where candidate draws are generated.

* `adapt` : type of adaptation phase.  Options are
    * `:all` : adapt proposal during all iterations.
    * `:burnin` : adapt proposal during burn-in iterations.
    * `:none` : no adaptation (Metropolis-within-Gibbs sampling with fixed proposal).
	

* `args...`: additional keyword arguments to be passed to the `AMWGVariate` constructor.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.AMWGVariate, ::Function; ::Bool=true)
```

#### Arguments

* `v`: current state of parameters to be simulated. When running the sampler in adaptive mode, the `v` argument in a successive call to the function will contain the `tune` field returned by the previous call.

* `adapt`: whether to adaptively update the proposal distribution.

## bhmc.jl

```@docs
Main.MCPhylo.BHMC
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.

* `traveltime`: length of time over which particle paths are simulated. It is recommended that supplied values be of the form ``(n + \frac{1}{2}) \pi``, where optimal choices of n ``\in \mathbb{Z}^+`` are expected to grow with the parameter space dimensionality. 

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.BHMCVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## bia.jl

```@docs
Main.MCPhylo.BIA
```

#### Arguments
* `params`: stochastic node(s) to be updated with the sampler.

* `args...`: additional keyword arguments to be passed to the `BIAVariate` constructor.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.BIAVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## bmc3.jl

```@docs
Main.MCPhylo.BMC3
```

#### Arguments 

* `params`: stochastic node(s) to be updated with the sampler.

* `k`: number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.SamplerVariate{Main.MCPhylo.BMC3Tune{Main.MCPhylo.BMC3Form}}, ::Function)
```

#### Arguments 

* `v`: current state of parameters to be simulated.

## bmg.jl

```@docs
Main.MCPhylo.BMG
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.

* `k`:  number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.SamplerVariate{Main.MCPhylo.BMGTune{Main.MCPhylo.BMGForm}}, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## dgs.jl

```@docs
Main.MCPhylo.DGS
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.DGSVariate, ::Function)
Main.MCPhylo.sample!(::Main.MCPhylo.DiscreteVariate, ::Vector{Float64})
```

#### Arguments

* `v`: current state of parameters to be simulated.

## hmc.jl

```@docs
Main.MCPhylo.HMC
```

#### Arguments 

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `epsilon`: step size.

* `L`: number of steps to take in the Leapfrog algorithm.

* `Sigma`: covariance matrix for the multivariate normal proposal distribution. The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated. If omitted, the identity matrix is assumed.

* `dtype` : differentiation for gradient calculations. Options are
    * `:central` : central differencing
    * `:forward` : forward differencing.
	
```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.HMCVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## mala.jl

```@docs
Main.MCPhylo.MALA
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `epsilon`: factor by which the drift and covariance matrix of the proposal distribution are scaled.

* `Sigma`: covariance matrix for the multivariate normal proposal distribution. The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated. If omitted, the identity matrix is assumed.

* `dtype` : differentiation for gradient calculations. Options are
    * `:central` : central differencing
    * `:forward` : forward differencing.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.MALAVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## miss.jl

```@docs
Main.MCPhylo.MISS
```

#### Arguments

* `params`: stochastic node(s) that contain missing values (`NaN`) to be updated with the sampler.

## nuts.jl

```@docs
Main.MCPhylo.NUTS
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `dtype` : differentiation for gradient calculations. Options are
    * `:central` : central differencing
    * `:forward` : forward differencing.
	
* `args...`: additional keyword arguments to be passed to the `NUTSVariate` constructor.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.NUTSVariate, ::Function; ::Bool=false)
```

#### Arguments

* `v`: current state of parameters to be simulated. When running the sampler in adaptive mode, the `v` argument in a successive call to the function will contain the `tune` field returned by the previous call.

* `adapt`: whether to adaptively update the `epsilon` step size parameter.

## rwm.jl

```@docs
Main.MCPhylo.RWM(::Main.MCPhylo.ElementOrVector{Symbol}, ::Main.MCPhylo.ElementOrVector{T}; args...) where {T<:Real}
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `scale`: scaling value or vector of the same length as the combined elements of nodes `params` for the `proposal` distribution. Values are relative to the unconstrained parameter space, where candidate draws are generated.

* `args...`: additional keyword arguments to be passed to the `RWMVariate` constructor.


```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.RWMVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

## slice.jl

```@docs
Main.MCPhylo.Slice
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.

* `width`: scaling value or vector of the same length as the combined elements of nodes `params`, defining initial widths of a hyperrectangle from which to simulate values.

* `F` : sampler type. Options are
    * `:Univariate` : sequential univariate sampling of parameters.
    * `:Multivariate` : joint multivariate sampling.

* `transform`: whether to sample parameters on the link-transformed scale (unconstrained parameter space). If `true`, then constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function, and `width` is interpreted as being relative to the unconstrained parameter space. Otherwise, sampling is relative to the untransformed space.

## slicesimplex.jl

```@docs
Main.MCPhylo.SliceSimplex
```

#### Arguments

* `params`: stochastic node(s) to be updated with the sampler.

* `args...`: additional keyword arguments to be passed to the `SliceSimplexVariate` constructor.

```@docs
Main.MCPhylo.sample!(::Main.MCPhylo.SliceSimplexVariate, ::Function)
```

#### Arguments

* `v`: current state of parameters to be simulated.

```@autodocs
Modules = [MCPhylo]
Pages   = ["sampler/rwmc.jl", "sampler.jl"]
Filter = 
```
