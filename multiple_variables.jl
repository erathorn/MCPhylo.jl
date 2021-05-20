using Revise
using Pkg
Pkg.activate(".")
using MCPhylo

data = Array{Float64,2}(undef, 15, 100)
μ1 = randn()
σ1 = rand(Exponential())
l = rand(LogNormal())
for i in 1:15
    μ = rand(Normal(μ1, σ1))
    σ = rand(Exponential(l))
    data[i,:] .= rand(Normal(μ, σ), 100)
end

my_data = Dict{Symbol, Any}(
  :df => data,
);


model = Model(
    df = Stochastic(2, (mu, sig) -> UnivariateDistribution[Normal(mu[i], sig[i]) for i in 1:15, j in 1:100], false, false),
    mu = Stochastic(1, (mu_hat, sig_hat)->Normal(mu_hat, sig_hat),true),
    sig = Stochastic(1, (l)->Exponential(l),true),
    l = Stochastic(()->LogNormal(),true),
    mu_hat = Stochastic(()->Normal(), true),
    sig_hat = Stochastic(()->Exponential(), true)
   );

inits = [ Dict{Symbol, Union{Any, Real}}(
    :df => my_data[:df],
    :mu => randn(15),
    :sig => rand(15),
    :mu_hat => randn(),
    :l => rand(),
    :sig_hat => rand(),
    ) for i in 1:2]

scheme = [
    Slice([:mu, :sig], 1.0),
    Slice([:mu_hat, :sig_hat, :l], 1.0)
]

setsamplers!(model, scheme)

sim = mcmc(model, my_data, inits, 1000, burnin=500, thin=10, chains=2)


# warning example
p = plot(sim; force=true, fuse=true, vars=["mu", "sig"], ymirror=true)
# example with the layout use
p = plot(sim; fuse_layout=(1,2), layout=(5,6), force=true, fuse=true, vars=["mu", "sig"])
p = plot(sim; fuse_layout=(2,1), layout=(5,6), force=true, fuse=true, vars=["mu", "sig"])
p[2][:size]
