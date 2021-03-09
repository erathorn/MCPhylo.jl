#################### MCMC Simulation Engine ####################


"""
    function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true, trees::Bool=false)

This function simulates additional draws from a model. `mc` is the results of
a previous call to the mcmc function. `iters` indicates the number of draws to
simulate, while `verbose` controls whether to print porgress statements to the
console. `trees` indicates if the states of the model
nodes describing tree structures should be stored as well.
"""
function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true, trees::Bool=false)
  thin = step(mc)
  last(mc) == div(mc.model.iter, thin) * thin ||
    throw(ArgumentError("chain is missing its last iteration"))

  mm = deepcopy(mc.model)
  mc2 = mcmc_master!(mm, mm.iter .+ (1:iters), last(mc), thin, mc.chains,
                     verbose, trees)
  if mc2.names != mc.names
    mc2 = mc2[:, mc.names, :]
  end

  ModelChains(vcat(mc, mc2), mc2.model)
end


"""
    mcmc(m::Model, inputs::Dict{Symbol},
              inits::Vector{V} where V<:Dict{Symbol},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true, trees::Bool=false)

Simulate MCMC draws from the model `m`. `inputs` is a dictionary storing the values
for input model nodes. Dictionary keys and values should be given for each input node.
The vector `inits` contains dictionaries with initial values for stochastic model nodes.
Dictionary keys and values should be given for each stochastic node. Consecutive
runs of the simulator will iterate through the vector’s dictionary elements.
The number of draws to generate for each simulation run is specified using the
`iter` argument. `burnin` specifies the number of initial draws to discard as a burn-in
sequence to allow for convergence. The thinning parameter `thin` is the step-size
between draws to store in the output. The parameter `chains` specifies the number
of simulation runs to perform. `verbose` indicates whether the sampler progress
should be printed at the console. `trees` indicates if the states of the model
nodes describing tree structures should be stored as well.
"""
function mcmc(m::Model, inputs::Dict{Symbol},
              inits::Vector{V} where V<:Dict{Symbol},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true, trees::Bool=false)
  iters > burnin ||
    throw(ArgumentError("burnin is greater than or equal to iters"))
  length(inits) >= chains ||
    throw(ArgumentError("fewer initial values than chains"))
  if trees == true
    nt = 0
    for n in values(m.nodes)
      if isa(n, TreeStochastic)
        nt += 1
      end
    end
    if nt == 0
      throw(ArgumentError("no tree obejct to write"))
    elseif nt > 1
      throw(ArgumentError("to many tree obejcts to write"))
    end
  end

  mm = deepcopy(m)
  setinputs!(mm, inputs)
  setinits!(mm, inits[1:chains])
  mm.burnin = burnin

  mcmc_master!(mm, 1:iters, burnin, thin, 1:chains, verbose, trees)
end


function mcmc_master!(m::Model, window::UnitRange{Int}, burnin::Integer,
                      thin::Integer, chains::AbstractArray{Int}, verbose::Bool, trees::Bool)
  states = m.states
  m.states = ModelState[]

  N = length(window)
  K = length(chains)

  frame = ChainProgressFrame(
    "MCMC Simulation of $N Iterations x $K Chain" * "s"^(K > 1), verbose
  )

  lsts = [
    Any[m, states[k], window, burnin, thin, ChainProgress(frame, k, N), trees]
    for k in chains
  ]
  results = pmap2(mcmc_worker!, lsts)

  sims  = Chains[results[k][1] for k in 1:K]
  model = results[1][2]
  model.states = ModelState[results[k][3] for k in sortperm(chains)]

  ModelChains(cat(sims..., dims=3), model)
end


function mcmc_worker!(args::Vector)
  m::Model, state::ModelState, window::UnitRange{Int}, burnin::Integer, thin::Integer, meter::ChainProgress, store_trees::Bool = args
  llname::AbstractString = "likelihood"
  treeind = 1
  simind = 1
  m.iter = first(window) - 1
  relist!(m, state.value)
  settune!(m, state.tune)

  pnames = vcat(names(m, true), llname)

  sim = Chains(last(window), length(pnames), start=burnin + thin, thin=thin,
               names=pnames)
  #sim = Chains(last(window), length(pnames), start=thin, thin=thin,
  #            names=pnames)
  treenode = :tn
  for i in m.nodes
    if isa(i[2], TreeStochastic)
      treenode = i[1]
      break
    end
  end

  reset!(meter)
  for i in window

    sample!(m)
    if i > burnin && (i - burnin) % thin == 0

      sim[i, :, 1] = unlist(m, true)

      if store_trees
       sim.trees[treeind, 1, 1] = newick(m[treenode].value)
       treeind +=1
     end
    end
    next!(meter)
  end

  mv = samparas(m)

  sim.moves[1] = mv
  (sim, m, ModelState(unlist(m), gettune(m)))
end

function track(i::Integer, burnin::Integer, thin::Integer,
               sim::Chains, m::Model, store_trees::Bool, treeind::Integer, treenode)
  if i > burnin && (i - burnin) % thin == 0
    sim[i, :, 1] = unlist(m, true)
    if store_trees
     sim.trees[treeind, 1, 1] = newick(m[treenode].value)
     treeind +=1
   end
  end
end
