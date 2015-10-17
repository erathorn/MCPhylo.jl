#################### Model-Based Posterior Statistics ####################

function dic(mc::ModelChains)
  m = mc.model
  nodekeys = keys(m, :output)
  idx = indexin(names(m, keys(m, :block)), mc.names)
  0 in idx && error("dic requires all sampled nodes to be monitored")

  xbar = map(i -> mean(mc.value[:,i,:]), idx)
  relist!(m, xbar)
  Dhat = -2.0 * mapreduce(key -> logpdf(m[key]), +, nodekeys)
  D = -2.0 * logpdf(mc, nodekeys)
  p = [mean(D) - Dhat, 0.5 * var(D)]

  ChainSummary([Dhat + 2.0 * p  p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(mc))
end

function logpdf(mc::ModelChains, nodekeys::Vector{Symbol})
  m = mc.model
  idx = indexin(names(m, keys(m, :block)), mc.names)
  0 in idx && error("logpdf requires all sampled nodes to be monitored")

  iters, p, chains = size(mc.value)
  values = Array(Float64, iters, 1, chains)
  frame = ChainProgressFrame(
    "MCMC Processing of $iters Iterations x $chains Chain" * "s"^(chains > 1),
    true
  )
  for k in 1:chains
    meter = ChainProgress(frame, k, iters)
    for i in 1:iters
      relist!(m, mc.value[i,idx,k])
      values[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nodekeys)
      next!(meter)
    end
    println()
  end

  values
end

function predict(mc::ModelChains, nodekey::Symbol)
  m = mc.model
  node = m[nodekey]

  outputs = keys(m, :output)
  nodekey in outputs ||
    error("predict is only defined for observed Stochastic nodes: ",
          join(map(string, outputs), ", "))

  nodenames = names(m, [nodekey])

  sources = intersect(node.sources, keys(m, :stochastic))
  idx = indexin(names(m, sources), mc.names)
  0 in idx && error("predict requires monitoring of nodes: ",
                      join(map(string, sources), ", "))

  iters, _, chains = size(mc.value)
  value = Array(Float64, iters, length(nodenames), chains)
  for k in 1:chains
    for i in 1:iters
      relist!(m, mc.value[i,idx,k], sources)
      value[i,:,k] = unlist(node, rand(node))
    end
  end

  Chains(value, start=start(mc.range), thin=step(mc.range), names = nodenames)
end
