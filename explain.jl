## Use Revise when editing the source code of the package
using Revise

using Pkg
Pkg.activate(".")

using MCPhylo

using Distributions


## simple function to extract the maximum point scale reduction factor
function max_psrf(sim::ModelChains)::Float64
    bi = 1 + size(sim)[1] ÷ 2
    gd = gelmandiag(sim[bi:end,:,:])
    gd_values::Array{Float64} = gd.value
    indices = isnan.(gd_values)
    gd_values[indices] .= -1
    psrf = maximum(gd_values[:,1], )
    psrf
end


## Generate Gold Standard

# store targets for future reference
target_mean = 0.0
target_sd = 1.0

gold_standard = Normal(target_mean, target_sd)
generated_data = rand(gold_standard, 1000)

## Setup Model
my_data = Dict{Symbol,Any}(
    :data => generated_data
)

model = Model(
    data = Stochastic(1, (me, sd)->Normal(me, sd), false),
    me = Stochastic(()->Normal(), true),
    sd = Stochastic(()->Exponential(), true)
)


initial_values = [Dict{Symbol,Any}(
    :data => generated_data,
    :me => rand(),
    :sd => rand()
) for i in 1:2]

# Select the sampler
scheme = [NUTS([:me, :sd])
            ]

setsamplers!(model, scheme)

## Run the Model
sim = mcmc(model, my_data, initial_values, 500, burnin=250,thin=10, chains=2)

psrf=max_psrf(sim)
while psrf > 1.1
      global sim, psrf
      @show psrf
      sim = mcmc(sim, 1000)
      psrf = max_psrf(sim)
end
