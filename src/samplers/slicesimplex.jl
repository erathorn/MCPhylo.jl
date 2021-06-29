#################### Slice Simplex Sampler ####################

#################### Types and Constructors ####################

mutable struct SliceSimplexTune <: SamplerTune
  logf::Union{Function, Missing}
  scale::Float64

  SliceSimplexTune() = new()

  SliceSimplexTune(x::Vector, logf::Union{Function, Missing}; scale::Real=1.0) =
    new(logf, scale)
end

SliceSimplexTune(x::Vector; args...) =
  SliceSimplexTune(x, missing; args...)

const SliceSimplexVariate = SamplerVariate{SliceSimplexTune}

function validate(v::SliceSimplexVariate)
  0 < v.tune.scale <= 1 || throw(ArgumentError("scale is not in (0, 1]"))
  validatesimplex(v)
end


#################### Sampler Constructor ####################
"""
    SliceSimplex(params::ElementOrVector{Symbol}; args...)

Construct a `Sampler` object for which SliceSimplex sampling is to be applied
separately to each of the supplied parameters. Parameters are assumed to be
continuous and constrained to a simplex.

Returns a `Sampler{SliceSimplexTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler.

* `args...`: additional keyword arguments to be passed to the `SliceSimplexVariate` constructor.
"""
function SliceSimplex(params::ElementOrVector{Symbol}; args...)
  params = asvec(params)
  samplerfx = function(model::Model, block::Integer)
    s = model.samplers[block]
    local node, x
    for key in params
      node = model[key]
      x = unlist(node)

      sim = let x = x, s=s, args = args
        function(inds::AbstractRange, logf::Function)
          v = SamplerVariate(x[inds], s, model.iter; args...)
          sample!(v, logf)
        end
      end

      logf = let x = x, key = key, model=model,node=node
        function(d::MultivariateDistribution, v::AbstractVector,
                      inds::AbstractRange)
          x[inds] = v
          relist!(model, x, key)
          logpdf(d, v) + logpdf(model, node.targets)
        end
      end

      SliceSimplex_sub!(node.distr, sim, logf)
    end

    nothing
  end
  Sampler(params, samplerfx, SliceSimplexTune())
end


function SliceSimplex_sub!(d::MultivariateDistribution, sim::Function,
                           logf::Function)
  inds = 1:length(d)
  sim(inds, v -> logf(d, v, inds))
end

function SliceSimplex_sub!(D::Array{MultivariateDistribution}, sim::Function,
                           logf::Function)
  inds = 0:0
  for i in 1:length(D)
    d = D[i]
    inds = last(inds) .+ (1:length(d))
    sim(inds, v -> logf(d, v, inds))
  end
end

function SliceSimplex_sub!(d, sim::Function, logf::Function)
  throw(ArgumentError("unsupported distribution structure $(typeof(d))"))
end


#################### Sampling Functions ####################

sample!(v::SliceSimplexVariate) = sample!(v, v.tune.logf)
"""
    sample!(v::SliceSimplexVariate, logf::Function)

Draw one sample from a target distribution using the SliceSimplex sampler.
Parameters are assumed to be continuous and constrained to a simplex.

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::SliceSimplexVariate, logf::Function)
  p0 = logf(v.value) + log(rand())

  d = Dirichlet(fill!(similar(v), 1))
  ct = 0
  vertices = makefirstsimplex(v, v.tune.scale)
  vb = vertices \ v
  xb = rand(d)
  x = vertices * xb
  while any(x .< 0.0) || any(x .> 1.0) || logf(x) < p0
    if ct > 10
      return v
    end
    vertices = shrinksimplex(vb, xb, v, x, vertices)
    vb = vertices \ v
    xb = rand(d)
    x = vertices * xb
    ct += 1
  end
  v[:] = x

  v
end


function makefirstsimplex(x::AbstractVector{Float64}, scale::Real)
  vertices = Matrix{Float64}(I, length(x), length(x))
  vertices[:, 2:end] += (1.0 - scale) * (vertices[:, 1] .- vertices[:, 2:end])
  vertices .+ x .- vertices * rand(Dirichlet(fill!(similar(x), 1)))
end


function shrinksimplex(bx::AbstractVector{Float64}, bc::AbstractVector{Float64},
                       cx::AbstractVector{Float64}, cc::AbstractVector{Float64},
                       vertices::AbstractMatrix{Float64})
  for i in findall(bc .< bx)
    inds = [1:(i - 1); (i + 1):size(vertices, 2)]
    vertices[:, inds] += bc[i] * (vertices[:, i] .- vertices[:, inds])
    bc = vertices \ cc
  end

  vertices
end
