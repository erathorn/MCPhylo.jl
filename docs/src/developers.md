
# For Developers

On this page you can finde more detailed, developer oriented documentation on internal functions
which are not necessarily exported to the main API.

```@docs
MCPhylo.setinits!
MCPhylo.setmonitor!
MCPhylo.update!
MCPhylo.Base.keys
MCPhylo.Base.show
MCPhylo.showall
```

## Simulation

```@docs
MCPhylo.gettune
MCPhylo.settune!
MCPhylo.gradlogpdf
MCPhylo.gradlogpdf!
MCPhylo.logpdf
MCPhylo.logpdf!
MCPhylo.sample!(::Main.MCPhylo.Model, ::Integer)
MCPhylo.unlist
MCPhylo.relist
MCPhylo.relist!
MCPhylo.update!
MCPhylo.sample!
MCPhylo.ABC_sample
MCPhylo.Sampler
```

```@docs
MCPhylo.getindex
MCPhylo.setindex!
MCPhylo.Base.cat
MCPhylo.Base.keys(::Main.MCPhylo.AbstractChains)
MCPhylo.Base.show(::IO, ::Main.MCPhylo.AbstractChains)
MCPhylo.Base.size
```
