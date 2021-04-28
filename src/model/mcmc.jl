#################### MCMC Simulation Engine ####################


"""
    mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true, trees::Bool=false)

This function simulates additional draws from a model.

* `mc` is the results of a previous call to the mcmc function.

* `iters` indicates the number of draws to simulate.

* `verbose` controls whether to print progress statements to the console.

* `trees` indicates if the states of the model nodes describing tree structures should be stored as well.
"""
function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true,
              trees::Bool=false, ASDSF::Bool=false, ASDSF_freq::Int64=100,
              ASDSF_min_splits::Float64=0.1)::ModelChains

  thin = step(mc)
  last(mc) == div(mc.model.iter, thin) * thin ||
    throw(ArgumentError("chain is missing its last iteration"))

  mm = deepcopy(mc.model)
  mc2 = mcmc_master!(mm, mm.iter .+ (1:iters), last(mc), thin, mc.chains,
                     verbose, trees, ASDSF, ASDSF_freq, ASDSF_min_splits)
  if mc2.names != mc.names
    mc2 = mc2[:, mc.names, :]
  end
  ModelChains(vcat(mc, mc2), mc2.model, cat(mc.stats, mc2.stats, dims=1), mc.stat_names)
end


"""
    mcmc(m::Model, inputs::Dict{Symbol},
              inits::Vector{V} where V<:Dict{Symbol},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true, trees::Bool=false)

Simulate MCMC draws from the model `m`.

* `inputs` is a dictionary storing the values for input model nodes. Dictionary keys and values should be given for each input node.

* `inits` contains dictionaries with initial values for stochastic model nodes. Dictionary keys and values should be given for each stochastic node. Consecutive runs of the simulator will iterate through the vector’s dictionary elements.

* `iter` Specifies the number of draws to generate for each simulation run

* `burnin` specifies the number of initial draws to discard as a burn-in sequence to allow for convergence.

* `thin` is the step-size between draws to store in the output.

* `chains` specifies the number of simulation runs to perform.

* `verbose` indicates whether the sampler progress should be printed at the console.

* `trees` indicates if the states of the model nodes describing tree structures should be stored as well.
"""
function mcmc(m::Model, inputs::Dict{Symbol},
              inits::Vector{V} where V<:Dict{Symbol},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true, trees::Bool=false,
              ASDSF::Bool=false, ASDSF_freq::Int64=100,
              ASDSF_min_splits::Float64=0.1)::ModelChains

  ASDSF && !trees &&
   throw(ArgumentError("ASDSF can not be calculated without trees"))
  iters > burnin ||
    throw(ArgumentError("burnin is greater than or equal to iters"))
  length(inits) >= chains ||
    throw(ArgumentError("fewer initial values than chains"))

  mm::Model = deepcopy(m)
  setinputs!(mm, inputs)
  setinits!(mm, inits[1:chains])
  mm.burnin = burnin
  mcmc_master!(mm, 1:iters, burnin, thin, 1:chains, verbose, trees, ASDSF,
               ASDSF_freq, ASDSF_min_splits)
end


function mcmc_master!(m::Model, window::UnitRange{Int}, burnin::Integer,
                      thin::Integer, chains::AbstractArray{Int}, verbose::Bool,
                      trees::Bool, ASDSF::Bool, ASDSF_freq::Int64,
                      ASDSF_min_splits::Float64)::ModelChains

  states::Vector{ModelState} = m.states
  m.states = ModelState[]

  N = length(window)
  K = length(chains)

  frame::ChainProgressFrame = ChainProgressFrame(
    "MCMC Simulation of $N Iterations x $K Chain" * "s"^(K > 1), verbose
  )

  lsts = [
    Any[m, states[k], window, burnin, thin, ChainProgress(frame, k, N), trees]
    for k in chains
  ]
  results::Vector{Tuple{Chains, Model, ModelState}}, stats::Array{Float64, 2}, statnames::Vector{AbstractString} = assign_mcmc_work(mcmc_worker!, lsts, ASDSF, ASDSF_freq, ASDSF_min_splits)

  sims::Array{Chains}  = Chains[results[k][1] for k in 1:K]
  model::Model = results[1][2]
  model.states = ModelState[results[k][3] for k in sortperm(chains)]
  ModelChains(cat(sims..., dims=3), model, stats, statnames)
end


function mcmc_worker!(args::Vector, ASDSF_step::Int64=100,
                      rc::Union{Nothing, RemoteChannel}=nothing
                      )::Tuple{Chains, Model, ModelState}
  m::Model, state::ModelState, window::UnitRange{Int}, burnin::Integer, thin::Integer, meter::ChainProgress, store_trees::Bool = args
  llname::AbstractString = "likelihood"
  treeind::Int64 = 1
  m.iter = first(window) - 1

  relist!(m, state.value)

  settune!(m, state.tune)
  pnames = vcat(names(m, true), llname)
  treenodes = Symbol[]
  for i in m.nodes
    if isa(i[2], TreeStochastic)
      push!(treenodes, i[1])
    end
  end

  sim = Chains(last(window), length(pnames), start=burnin + thin, thin=thin,
               names=pnames, ntrees=length(treenodes), tree_names=treenodes)

  reset!(meter)
  for i in window

    sample!(m)
    if i > burnin
      if (i - burnin) % thin == 0
        sim[i, :, 1] = unlist(m, true)
        if store_trees
          for (ind, tree_node) in enumerate(treenodes)
            sim.trees[treeind, ind, 1] = newick(m[tree_node].value)
          end # for
          treeind +=1
        end # if
      end # if
      if !isnothing(rc) && (i - burnin) % ASDSF_step == 0
        trees::Vector{AbstractString} = []
        for (ind, tree_node) in enumerate(treenodes)
          push!(trees, newick(m[tree_node].value))
        end # for
        put!(rc, trees)
      end # if
    end # if
    next!(meter)

  end # for
  mv = samparas(m)
  sim.moves[1] = mv
  (sim, m, ModelState(unlist(m), gettune(m)))
end # mcmc_worker!


function track(i::Integer, burnin::Integer, thin::Integer,
               sim::Chains, m::Model, store_trees::Bool, treeind::Integer, treenode)
  if i > burnin && (i - burnin) % thin == 0
    sim[i, :, 1] = unlist(m, true)
    if store_trees
     sim.trees[treeind, 1, 1] = newick(m[treenode].value)
     treeind +=1
   end # if
 end # if
end # track
