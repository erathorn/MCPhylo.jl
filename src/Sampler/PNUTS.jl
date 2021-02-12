#################### Phylogenetic No-U-Turn Sampler ####################

#################### Types and Constructors ####################

mutable struct PNUTSTune <: SamplerTune
  logfgrad::Union{Function, Missing}
  adapt::Bool
  alpha::Float64
  epsilon::Float64
  epsilonbar::Float64
  gamma::Float64
  Hbar::Float64
  kappa::Float64
  m::Int
  mu::Float64
  nalpha::Int
  t0::Float64
  delta::Float64
  target::Float64
  moves::Int
  tree_depth::Int
  targetNNI::Int
  nni::Int


  PNUTSTune() = new()

  function PNUTSTune(x::Vector{T}, epsilon::Float64, logfgrad::Union{Function, Missing};
                    target::Real=0.6, tree_depth::Int=10, targetNNI::Int=5) where T<:Node
    new(logfgrad, false, 0.0, epsilon, 1.0, 0.05, 0.0, 0.75, 0, NaN, 0, 10.0,0.003,
        target,0, tree_depth,targetNNI,0)
  end
end

PNUTSTune(x::Vector{T}, logfgrad::Function, ::NullFunction, delta::Float64=0.003; args...) where T<:Node =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)

PNUTSTune(x::Vector{T}, logfgrad::Function, delta::Float64; args...) where T<:Node =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)

PNUTSTune(x::Vector; epsilon::Real, args...) =
    NUTSTune(x, epsilon, missing, args...)

const PNUTSVariate = SamplerVariate{PNUTSTune}


#################### Sampler Constructor ####################

function PNUTS(params::ElementOrVector{Symbol}; dtype::Symbol=:forward,args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    f = (x, sz, ll, gr) -> mlogpdfgrad!(block, x, sz, ll, gr)
    v = SamplerVariate(block, f, NullFunction(); args...)

    sample!(v::PNUTSVariate, f, adapt=model.iter <= model.burnin)

    relist(block, v)
  end
  Sampler(params, samplerfx, PNUTSTune())
end


#################### Sampling Functions ####################

function mlogpdfgrad!(block::SamplingBlock, x::T, sz::Int64, ll::Bool=false, gr::Bool=false)::Tuple{Float64, Vector{Float64}}  where T<:Node
  grad = Vector{Float64}(undef, sz)
  lp = zero(Float64)

  if gr
    lp, grad = gradlogpdf!(block, x)::Tuple{Float64, Vector{Float64}}
  elseif ll
    lp = logpdf!(block, x)::Float64
  end
  lp, grad
end
sample!(v::PNUTSVariate; args...) = sample!(v, v.tune.logfgrad; args...)

function sample!(v::PNUTSVariate, logfgrad::Function; adapt::Bool=false)

  tune = v.tune
  setadapt!(v, adapt)
  if tune.adapt
     tune.m += 1

     nuts_sub!(v, tune.epsilon, logfgrad)
     HTACC = (tune.target - tune.alpha / tune.nalpha)
     avg = (tune.nni+1)/tune.nalpha
     HTNNI = (1/tune.targetNNI - 1/avg)

     HT = (HTACC + HTNNI) / 2

     p = 1.0 / (tune.m + tune.t0)
     tune.Hbar = (1.0 - p) * tune.Hbar +
                 p * HT
     tune.epsilon = exp(tune.mu - sqrt(tune.m) * tune.Hbar / tune.gamma)
     p = tune.m^-tune.kappa
     tune.epsilonbar = exp(p * log(tune.epsilon) +
                           (1.0 - p) * log(tune.epsilonbar))
  else
     if (tune.m > 0) tune.epsilon = tune.epsilonbar end
     nuts_sub!(v, tune.epsilon, logfgrad)
  end
  v
end


function setadapt!(v::PNUTSVariate, adapt::Bool)
  tune = v.tune
  if adapt && !tune.adapt
    tune.m = 0
    tune.mu = log(10.0 * tune.epsilon)
  end
  tune.adapt = adapt
  v
end



function nuts_sub!(v::P, epsilon::Float64, logfgrad::Function)::P where P<:PNUTSVariate

  mt::T where T<:Node = v.value[1]
  nl = size(mt)[1]-1
  delta = v.tune.delta

  r::Vector{Float64} = randn(nl)

  x, r, logf, grad, nni = refraction(mt, r, 1, zeros(nl), epsilon, logfgrad, delta, nl)

  lu::Float64 = log(rand())

  logp0 = logf - 0.5 * dot(r)
  logu0 = logp0 + lu
  rminus = rplus = r
  gradminus = gradplus = grad
  xminus = xplus = x

  j = 0
  n = 1
  s = true

  while s && j < v.tune.tree_depth

    pm =2 * (rand() > 0.5) - 1

    if pm == -1

      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha, nni1, lpp, nnip = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0, delta, nl,lu)

    else

      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha, nni1, lpp, nnip=
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0, delta, nl, lu)

    end#if pm

    if sprime && rand() < nprime / n
        v.value[1]= xprime
    end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus, gradminus, gradplus, epsilon, logfgrad, delta, nl, j)
    v.tune.alpha, v.tune.nalpha = alpha, nalpha
    nni += nni1
    v.tune.nni = nnip
  end
  v.tune.moves += nni
  v
end


function refraction(v::T, r::V, pm::I,
                    grad::V, epsilon::F, logfgrad::Function,
                    delta::F, sz::I)::Tuple{T,V, F, V, I}  where {T<:Node, V<:Vector{Float64}, F<:Real, I<:Integer}

    ref_r::V = pm*r

    blenvec::V = get_branchlength_vector(v)
    fac::V = scale_fac.(blenvec, delta)


    ref_r = @. ref_r + (epsilon * 0.5) * (grad * fac)
    tmpB::V = @. blenvec + (epsilon * ref_r)

    nni::Int64 = 0

    if minimum(tmpB) <= 0.0
        v, tmpB, ref_r, nni = ref_NNI(v, tmpB, ref_r, epsilon, blenvec, delta, logfgrad, sz)

    end

    blenvec = molifier.(tmpB, delta)

    set_branchlength_vector!(v, blenvec)


    logf, grad = logfgrad(v, sz, true, true)

    fac = scale_fac.(blenvec, delta)
    ref_r = @. ref_r + (epsilon * 0.5) * (grad * fac)

    return v, ref_r, logf, grad, nni
end


function ref_NNI(v::T, tmpB::Vector{Float64}, r::V, epsilon::F, blv::V,
                 delta::F, logfgrad::Function, sz::I)::Tuple{T, V, V, I}  where {T<:Node, V<:Vector{Float64}, I<:Integer, F<:Real}

  intext = internal_external(v)
  t = zero(F)
  nni = zero(I)

  while minimum(tmpB)<=0.0

     timelist = tmpB./abs.(r)
     ref_index::Int64 = argmin(timelist)
     temp::Float64 = epsilon-t+timelist[ref_index]
     blv = @. blv + (temp * r)

     r[ref_index] *= -1.0

     if Bool(intext[ref_index])

       blv1::Vector{Float64} = molifier.(blv, delta)
       set_branchlength_vector!(v, blv1)

       # use thread parallelism
       res_before = @spawn logfgrad(v, sz, true, false) # still with molified branch length

       v_copy = deepcopy(v)
       tmp_NNI_made = NNI!(v_copy, ref_index)

       # fetch the results from the parallel part
       U_before_nni, _ = fetch(res_before)
       U_before_nni *= -1

       if Bool(tmp_NNI_made)
            U_after_nni, _ = logfgrad(v_copy, sz, true, false)
            U_after_nni *= -1
            delta_U::Float64 = 2.0*(U_after_nni - U_before_nni)
            my_v::Float64 = r[ref_index]^2
            if my_v > delta_U
              nni += tmp_NNI_made
              r[ref_index] = sqrt(my_v - delta_U)
              v = v_copy
            end # if my_v
       end #if NNI
      end #non leave
      t = epsilon+timelist[ref_index]

      tmpB = @. blv + (epsilon-t) * r

  end #while

  v, tmpB, r, nni
end



function buildtree(x::T, r::V,
                   grad::V, pm::I, j::I,
                   epsilon::F, logfgrad::Function, logp0::F, logu0::F,
                   delta::F, sz::I, lu::F)::Tuple{T, V, V, T, V, V, T, I, Bool, F, F, I, F, I}  where {T<:Node, V<:Vector{Float64}, I<:Integer, F<:Real}


  if j == 0

    xprime, rprime, logfprime, gradprime, nni = refraction(x, r, pm, grad, epsilon,
                                          logfgrad, delta, sz)

    logpprime = logfprime - 0.5 * dot(rprime)

    nprime = Int(logu0 < logpprime)

    sprime = logu0 < logpprime + 1000.0
    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    nalphaprime = 1
    nnip = nni
  else
    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime, nni, logpprime, nnip = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0, delta,  sz, lu)
    if sprime

      if pm == -1

        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 , nni, logpprimem, nnip2 = buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta,  sz, lu)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2, nni, logpprime, nnip2 = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta, sz,lu)
      end # if pm

      if rand() < nprime2 / (nprime + nprime2)
        xprime = xprime2
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus, gradminus, gradplus, epsilon, logfgrad, delta, sz, j)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
      nnip += nnip2
    end #if sprime
  end #if j

  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime, nni, logpprime, nnip
end


function nouturn(xminus::T, xplus::T,
                rminus::V, rplus::V, gradminus::V,gradplus::V,
                epsilon::R, logfgrad::Function, delta::R, sz::I, j::I)::Bool  where {T<:Node, V<:Vector{Float64}, R <: Real, I <: Int}

        curr_l, curr_h = BHV_bounds(xminus, xplus)

        # use thread parallelism to calculuate both directions at once
        res_minus = Base.Threads.@spawn refraction(deepcopy(xminus), deepcopy(rminus), -1, gradminus, epsilon, logfgrad, delta, sz)
        xplus_bar,_,_,_,_ = refraction(deepcopy(xplus), deepcopy(rplus),1, gradplus, epsilon, logfgrad, delta,sz)

        # fetch the results
        xminus_bar,_,_,_,_ = fetch(res_minus)

        curr_t_l, curr_t_h = BHV_bounds(xminus_bar, xplus_bar)
        return curr_h < curr_t_l
end




#################### Auxilliary Functions ####################

function nutsepsilon(x::Node, logfgrad::Function, delta::Float64)

  x0 = deepcopy(x)
  n = size(x)[1] - 1

  _, r0, logf0, grad0,_ = refraction(x0, randn(n), 1, zeros(n), 0.0, logfgrad, delta, n)

  x0 = deepcopy(x)
  epsilon = 1.0
  _, rprime, logfprime, gradprime,_ = refraction(x0, r0, 1, grad0,  epsilon, logfgrad, delta, n)

  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  pm = 2 * (prob > 0.5) - 1
  while prob^pm > 0.5^pm
    epsilon *= 2.0^pm
    x0 = deepcopy(x)
    _, rprime, logfprime, _ ,_ = refraction(x0, r0, 1, grad0, epsilon, logfgrad, delta, n)

    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  end

  epsilon
end

@inline function scale_fac(x::T, delta::T)::T where T<:Real
  x < delta ? x/delta : 1.0
end

@inline function molifier(x::T, delta::T)::T where T<:Real
    x >= delta ? x : (x^2+delta^2)/(2.0*delta)
end # function
