#################### Model Simulation ####################

function gettune(m::Model, block::Integer)
  block == 0 && return gettune(m)
  m.samplers[block].tune
end

function gettune(m::Model)
  Any[gettune(m, i) for i in 1:length(m.samplers)]
end

function settune!(m::Model, tune, block::Integer)
  block == 0 && return settune!(m, tune)
  m.samplers[block].tune = tune
end

function settune!(m::Model, tune::Vector{Any})
  nsamplers = length(m.samplers)
  ntune = length(tune)
  nsamplers == ntune ||
    throw(DimensionMismatch(
      "tried to assign $ntune tune elements to $nsamplers samplers"
    ))

  for i in 1:nsamplers
    settune!(m, tune[i], i)
  end
end


function gradlogpdf(m::Model, block::Integer=0, transform::Bool=false;
                    dtype::Symbol=:forward)
  x0 = unlist(m, block, transform)
  value = gradlogpdf!(m, x0, block, transform, dtype=dtype)
  relist!(m, x0, block, transform)
  value
end


function gradlogpdf!(m::Model, x::AbstractVector{T}, block::Integer=0,
                      transform::Bool=false; dtype::Symbol=:forward) where {T<:Real}
  f = x -> logpdf!(m, x, block, transform)
  if dtype == :Zygote
    Zygote.gradient(f, x)
  else
    FiniteDiff.finite_difference_gradient(f, x)
  end

end


function logpdf(m::Model, block::Integer=0, transform::Bool=false)
  params = keys(m, :block, block)
  targets = keys(m, :target, block)
  logpdf(m, params, transform) + logpdf(m, setdiff(targets, params))
end

function logpdf(m::Model, nodekeys::Vector{Symbol}, transform::Bool=false)
  lp = Base.Threads.Atomic{Float64}(0.0)
  Threads.@threads for key in nodekeys
    Base.Threads.atomic_add!(lp, logpdf(m[key], transform))
    #lp += logpdf(m[key], transform)
    isfinite(lp[]) || break
  end

  m.likelihood=lp[]
  m.likelihood
end

function rand(m::Model, nodekeys::Vector{Symbol}, x::Int64)
  rand(m[nodekeys[1]], x)
end



function logpdf(m::Model, x::AbstractArray{T}, block::Integer=0,
                transform::Bool=false) where {T<:Real}
  x0 = unlist(m, block)
  lp = logpdf!(m, x, block, transform)
  relist!(m, x0, block)
  lp
end


function logpdf!(m::Model, x::N, block::Integer=0,
                  transform::Bool=false) where N<:GeneralNode
  params = keys(m, :block, block)
  targets = keys(m, :target, block)
  m[params] = relist(m, x, params, transform)
  lp = logpdf(m, setdiff(params, targets), transform)
  for key in targets
    isfinite(lp) || break
    node = m[key]
    update!(node, m)
    lp += key in params ? logpdf(node, transform) : logpdf(node)
  end
  lp
end

function gradlogpdf!(m::Model, x::AbstractArray{T}, block::Integer=0,transform::Bool=false) where T<:GeneralNode
  gradlogpdf!(m, x, block, transform)
end


function gradlogpdf(m::Model, targets::Array{Symbol, 1})::Tuple{Float64, Array{Float64}}
  vp = 0.0
  gradp = Array[]
  for key in targets
    node = m[key]
    update!(node, m)
    v, grad = gradlogpdf(node)
    vp += v
    push!(gradp, grad)
  end
  m.likelihood = vp
  vp, .+(gradp...)
end

function gradlogpdf!(m::Model, x::N, block::Integer=0,transform::Bool=false)::Tuple{Float64, Vector{Float64}} where N<:GeneralNode
  params = keys(m, :block, block)
  targets = keys(m, :target, block)
  m[params] = relist(m, x, params, transform)


  # use thread parallelism
  # likelihood
  lik_res = @spawn gradlogpdf(m, targets)
  
  # prior
  vp, gradp =  gradlogpdf(m[params[1]], x)

  # get results from threads
  v, grad = fetch(lik_res)
  vp+v, gradp.+grad
end

function logpdf!(m::Model, x::AbstractArray{T}, block::Integer=0,
                  transform::Bool=false) where {T<:Real}
  params = keys(m, :block, block)
  targets = keys(m, :target, block)
  m[params] = relist(m, x, params, transform)
  lp = logpdf(m, setdiff(params, targets), transform)
  for key in targets
    isfinite(lp) || break
    node = m[key]
    update!(node, m)
    lp += key in params ? logpdf(node, transform) : logpdf(node)
  end
  lp
end

function rand!(m::Model, x::Int64, block::Integer=0)
  params = keys(m, :block, block)
  res = rand(m, params, x)
  res
end


function sample!(m::Model, block::Integer=0)
  m.iter += 1
  isoneblock = block != 0
  blocks = isoneblock ? block : 1:length(m.samplers)
  for b in blocks
    sampler = m.samplers[b]
    value = sampler.eval(m::Model, b::Int)
    if value != nothing
      m[sampler.params] = value
      update!(m, b)
    end
  end
  m.iter -= isoneblock
  m.likelihood = logpdf(m)
  m
end


function final_likelihood(model::Model)::Float64
  logpdf(model[keys_output(model)[1]])
end

function unlist(m::Model, block::Integer=0, transform::Bool=false)
  unlist(m, keys(m, :block, block), transform)
end

function samparas(m::Model)
  moves = 0
  for i in m.samplers
    if typeof(i.tune) == PNUTSTune
      moves = i.tune.moves
    end
  end
 return moves
end

function unlist(m::Model, monitoronly::Bool)
  f = function(key)
    node = m[key]
    lvalue = isa(node, AbstractTreeStochastic) ? unlist_tree(node) : unlist(node)
    monitoronly ? lvalue[node.monitor] : lvalue
  end
  r = vcat(map(f, keys(m, :dependent))..., m.likelihood)
  r
end

function unlist(m::Model, nodekeys::Vector{Symbol}, transform::Bool=false)
  vcat(map(key -> unlist(m[key], transform), nodekeys)...)
end


function relist(m::Model, x::AbstractArray{T}, block::Integer=0,
                transform::Bool=false) where {T<:Real}
  relist(m, x, keys(m, :block, block), transform)
end

function relist(m::Model, x::AbstractArray{T}, block::Integer=0,
               transform::Bool=false) where {T<:GeneralNode}

 relist(m, x, keys(m, :block, block), transform)
end


function relist(m::Model, x::AbstractArray{T},
                nodekeys::Vector{Symbol}, transform::Bool=false) where {T<:Any}
  values = Dict{Symbol,Union{Any, Real}}()

  N = length(x)
  offset = 0
  for key in nodekeys

    value, n = relistlength(m[key], view(x, (offset + 1):N), transform)

    values[key] = value
    offset += n
  end
  offset == length(x) ||
    throw(ArgumentError("incompatible number of values to put in nodes"))
  values
end



function relist(m::Model, x::N, nodekeys::Vector{Symbol}, transform::Bool=false) where N<:GeneralNode
  values = Dict{Symbol,Any}()
  offset = 0
  for key in nodekeys
    value, n = relistlength(m[key], x, transform)
    values[key] = value
    offset += n
  end

  values
end



function relist!(m::Model, x::AbstractArray{T}, block::Integer=0,
                 transform::Bool=false) where {T<:Any}
  nodekeys = keys(m, :block, block)
  values = relist(m, x, nodekeys, transform)
  for key in nodekeys
    assign!(m, key, values[key])
  end
  update!(m, block)
end

function assign!(m::Model, key::Symbol, value::T) where T <: Real
  m[key].value = value
end

function assign!(m::Model, key::Symbol, value::T) where T <:GeneralNode
  m[key].value = value
end

function assign!(m::Model, key::Symbol, value::T) where T <: Array
  if isa(m[key], TreeVariate)
    @assert (2,) == size(value)
    d = Array{Float64, 1}(undef, 2)
  else
    @assert size(m[key].value) == size(value)
    d = similar(m[key].value)
  end
  for (ind, elem) in enumerate(value)
    d[ind] = elem
  end
  m[key].value = d
end



function relist!(m::Model, x::AbstractArray{T}, nodekey::Symbol,
                  transform::Bool=false) where {T<:Real}
  node = m[nodekey]
  m[nodekey] = relist(node, x, transform)
  update!(m, node.targets)
end


function update!(m::Model, block::Integer=0)
  nodekeys = block == 0 ? keys(m, :dependent) : m.samplers[block].targets
  update!(m, nodekeys)
end

function update!(m::Model, nodekeys::Vector{Symbol})
  for key in nodekeys
    update!(m[key], m)
  end
  m
end
