using Revise
using Pkg
Pkg.activate(".")
using MCPhylo


using Serialization


data = rand(Normal(0,1), 5000)

my_data=Dict(:data=>data)

model = Model(
     data = Stochastic(1, (μ, σ) -> Normal(μ, σ), false),
        μ = Stochastic(()->Normal(),true),
        σ = Stochastic(()->Exponential(1), true)
)

inits = [Dict(:data => data,
             :μ => randn(),
             :σ => rand()),
        Dict(:data => data,
            :μ => randn(),
            :σ => rand())]

samplers = [NUTS(:μ),
            Slice(:σ, 0.1)]

setsamplers!(model, samplers)

p = SimulationParameters()
p.asdsf = true

sim = mcmc(model, my_data, inits, 5000, burnin=500, thin=5, chains=2, params=p)

plot(sim, [:autocor, :contour, :bar, :mean, :density, :trace, :mixeddensity], fmt=:pdf, nrow=2, ncol=1)

gelmandiag(sim)

describe(sim)


# trees = ParseNewick("Drav_mytrees_1.nwk")
# trees = ParseNewick("gs_tree_JC_20-60-letters.nwk")

# t1 = deepcopy(trees[1])
# t2 = deepcopy(trees[2])
t1 = ParseNewick("(((((A:1,B:1):0.88,(C:1,D:1):1)J:0.47,E:1):0.73,F:1):0.83,G:1);")
t2 = ParseNewick("(((((C:0.2,D:1):0.5,A:1):0.15,E:1):0.87,(B:1,G:1):0.42):0.7,F:1);")

geo = geodesic(t1, t2; verbose=true)
println("Geodesic distance of the two trees is: $geo")

MCPhylo.ParseNexus("Example.nex")