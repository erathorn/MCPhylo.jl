mutable struct BranchSliceTune <: SamplerTune
    logf::Union{Function, Missing}
    width::Float64

    BranchSliceTune() = new()

    BranchSliceTune(width::Float64, func::Union{Missing,Function}) = new(func, width)
    function BranchSliceTune(x::Vector, width::Float64, func::Union{Function, Missing})
      new(func, width)
    end

end # mutable struct

BranchSliceTune(x::Vector,width::Float64 ; args...) =
  BranchSliceTune(x, width, missing; args...)


const BranchSliceVariate = SamplerVariate{BranchSliceTune}


 function BranchSlice(params::ElementOrVector{Symbol},
                      width::ElementOrVector{T},
                      transform::Bool=false) where {T<:Real}

      samplerfx = function(model::Model, block::Integer)

          block = SamplingBlock(model, block, transform)
          v = SamplerVariate(block, width)
          sample!(v, x -> logpdf!(block, x))
          relist(block, v)
      end # function samplerfx
      Sampler(params, samplerfx, BranchSliceTune())
  end

sample!(v::BranchSliceVariate) = sample!(v, v.tune.logf)

function sample!(v::BranchSliceVariate, logf::Function)
  p0 = logf(v) + log(rand())
  blens = get_branchlength_vector(v.value[1])

  n = length(blens)
  lower = blens - v.tune.width .* rand(n)
  upper = lower .+ v.tune.width

  x = v.tune.width .* rand(n) + lower
  set_branchlength_vector!(v.value[1], x)
  while logf(v) < p0
    for i in 1:n
      value = x[i]
      if value < blens[i]
        lower[i] = value
      else
        upper[i] = value
      end
      x[i] = rand(Uniform(lower[i], upper[i]))
    end
    set_branchlength_vector!(v.value[1], x)
  end
  #v[:] = x

  v
end
