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


  PNUTSTune() = new()

  function PNUTSTune(x::Vector{T}, epsilon::Float64, logfgrad::Union{Function, Missing};
                    target::Real=0.6) where T<:Node
    new(logfgrad, false, 0.0, epsilon, 1.0, 0.05, 0.0, 0.75, 0, NaN, 0, 10.0,0.003,
        target,0)
  end
end

PNUTSTune(x::Vector{T}, logfgrad::Function, ::NullFunction, delta::Float64=0.003; args...) where T<:Node =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)

PNUTSTune(x::Vector{T}, logfgrad::Function, delta::Float64; args...) where T<:Node =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)


  NUTSTune(x::Vector{T}, epsilon::Float64; args...) where T<:Node =
    NUTSTune(x, epsilon, missing; args...)

const PNUTSVariate = SamplerVariate{PNUTSTune}


#################### Sampler Constructor ####################

function PNUTS(params::ElementOrVector{Symbol}; dtype::Symbol=:forward, args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    f = (x, sz, ll, gr) -> mlogpdfgrad!(block, x, sz, ll, gr)
    v = SamplerVariate(block, f, NullFunction(); args...)

    sample!(v, f, adapt=model.iter <= model.burnin)
    relist(block, v)
  end
  Sampler(params, samplerfx, PNUTSTune())
end


#################### Sampling Functions ####################

function mlogpdfgrad!(block::SamplingBlock, x::T, sz::Int64, ll::Bool=false, gr::Bool=false)  where T<:Node
  grad = Vector{Float64}(undef, sz)
  lp = 0.0
  if gr
    lp, grad = gradf!(block, x, :provided)
  elseif ll
    lp = logpdf!(block, x)
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
    p = 1.0 / (tune.m + tune.t0)
    tune.Hbar = (1.0 - p) * tune.Hbar +
                p * (tune.target - tune.alpha / tune.nalpha)
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



function nuts_sub!(v::PNUTSVariate, epsilon::Float64, logfgrad::Function)

  mt = v.value[1]
  nl = size(mt)[1]-1
  delta = v.tune.delta

  #x = deepcopy(mt)

  #currU, _ = logfgrad(mt, nl, true, false)
  r = randn(nl)
  g = zeros(nl)

  lor = bitrand(nl)

  x, r , logf, grad, nni = refraction(mt, r, 1.0, g, epsilon, logfgrad, delta, lor, nl)

  #probH = logf-0.5*dot(r)
  #ratio = currH - probH
  lu = log(rand())

  logp0 = logf - 0.5 * dot(r)
  logu0 = logp0 + lu
  rminus = rplus = r
  gradminus = gradplus = grad
  #if ratio > min(0, log(rand()))
  xminus = xplus = x
  #else
    #logp0 = currU - 0.5 * dot(org_r)
    #logu0 = logp0 + log(rand())
  #  xminus = xplus = mt
  #  #rminus1 = rplus1 = org_r
    #gradminus1 = gradplus1 = grad
  #  nni = 0
  #end

  j = 0
  n = 1
  s = true

  while s


    pm =Float64( 2 * (rand() > 0.5) - 1)

    if pm == -1.0

      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha, nni1, lpp = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0, delta, lor, nl,lu)

    else

      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha, nni1, lpp =
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0, delta, lor, nl, lu)

    end#if pm
    #ratio = currH-lpp
    if sprime && rand() < nprime / n #&& rand() < exp(lpp-currH)
        v.value[1]= xprime
    end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus, gradminus, gradplus, epsilon, logfgrad, delta, nl, j,lor)
    v.tune.alpha, v.tune.nalpha = alpha, nalpha
    nni += nni1
  end
  v.tune.moves += nni
  v
end

@inline function scale_fac(x::T, delta::T) where T<:Float64
  x < delta ? x/delta : 1.0
end

function refraction(v::T, r::Vector{Float64}, pm::Float64,
                    grad::Vector{Float64}, epsilon::Float64, logfgrad::Function,
                    delta::Float64, lor::BitArray, sz::Int64)  where T<:Node

    v1 = deepcopy(v)
    ref_grad = grad
    ref_r = r
    if pm < 0
      mlor = .!lor
    else
      mlor = lor
    end

    blenvec = get_branchlength_vector(v1)
    fac = scale_fac.(blenvec, delta)

    ref_r = @. ref_r - (epsilon * 0.5) * (ref_grad * fac)

    tmpB = @. blenvec + (epsilon * ref_r)

    nni = 0

    if minimum(tmpB) <= 0
        v1, tmpB, ref_r, nni = ref_NNI(v1, tmpB, ref_r, epsilon, blenvec, delta, logfgrad, mlor, sz)

    end
    #println(tmpB)
    set_branchlength_vector!(v1, molifier.(tmpB, delta))
    #set_branchlength_vector!(v1, tmpB)
    logf, grad = logfgrad(v1, sz, true, true)
    #if pm > 0
    #  ref_grad = -grad
    #else
    #  ref_grad = grad
    #end
    get_branchlength_vector(v1, blenvec)
    fac = scale_fac.(blenvec, delta)
    ref_r = @. ref_r - (epsilon * 0.5) * (ref_grad * fac)
    #r = @. r - (epsilon * 0.5) * ref_grad# * fac)
    return v1, ref_r, logf, grad, nni
end


"""
    molifier(x::Float64, delta::Float64)::Float64

documentation
"""
@inline function molifier(x::Float64, delta::Float64)::Float64
    x >= delta ? x : (x^2+delta^2)/(2.0*delta)
end # function



function ref_NNI(v::T, tmpB::Vector{Float64}, r::Vector{Float64}, epsilon::Float64, blv::Vector{Float64},
                 delta::Float64, logfgrad::Function, lor::BitArray, sz::Int64)  where T<:Node

  intext = internal_external(v)
  t = 0.0
  nni = 0
  ## r = probM
  while minimum(tmpB)<=0.0

     timelist = tmpB./abs.(r)
     ref_index = argmin(timelist)

     temp = epsilon-t+timelist[ref_index]
     blv = @. blv + (temp * r)

     r[ref_index] *= -1.0

     if intext[ref_index] == 1

       blv1 =molifier.(blv, delta)
       set_branchlength_vector!(v, blv1)
       U_before_nni, _ = logfgrad(v, sz, true, false) # still with molified branch length

       U_before_nni *= -1
       v_copy = deepcopy(v)
       tmp_NNI_made = NNI!(v_copy, ref_index, lor[ref_index])
       nni += NNI!(v, ref_index, lor[ref_index])
       if tmp_NNI_made != 0

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



function buildtree(x::T, r::Vector{Float64},
                   grad::Vector{Float64}, pm::Float64, j::Integer,
                   epsilon::Float64, logfgrad::Function, logp0::Real, logu0::Real,
                   delta::Float64, lor::BitArray, sz::Int64, lu::Float64)  where T<:Node


  if j == 0
    #x1 = deepcopy(x)
    xprime, rprime, logfprime, gradprime, nni = refraction(x, pm*r, pm, grad, epsilon,
                                          logfgrad, delta, lor, sz)

    logpprime = logfprime - 0.5 * dot(rprime)

    nprime = Int(logu0 < logpprime)
    #nprime = Int(lu <= exp(logpprime))
    sprime = logu0 < logpprime + 1000.0

    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    nalphaprime = 1

  else
    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime, nni, logpprime = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0, delta, lor, sz, lu)
    if sprime
      if pm == -1.0

        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 , nni, logpprime= buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta, lor, sz, lu)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2, nni, logpprime = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta, lor, sz,lu)
      end # if pm

      if rand() < nprime2 / (nprime + nprime2)
        xprime = xprime2
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus, gradminus, gradplus, epsilon, logfgrad, delta, sz, j, lor)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
    end #if sprime
  end #if j

  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime, nni, logpprime
end


function nouturn(xminus::T, xplus::T,
                rminus::Vector{Float64}, rplus::Vector{Float64}, gradminus::Vector{Float64},gradplus::Vector{Float64},
                epsilon::Float64, logfgrad::Function, delta::Float64, sz::Int64, j::Int64, lor::BitArray)  where T<:Node

        curr_l, curr_h = BHV_lower(xminus, xplus)
        xminus_bar,_,_,_,_ = refraction(deepcopy(xminus), deepcopy(rminus), -1.0, gradminus, epsilon, logfgrad, delta, lor, sz)
        xplus_bar,_,_,_,_ = refraction(deepcopy(xplus), deepcopy(rplus),1.0, gradplus, epsilon, logfgrad, delta, lor, sz)
        curr_t_l, curr_t_h = BHV_lower(xminus_bar, xplus_bar)
        return curr_h < curr_t_l
          #return true
        #elseif curr_l > curr_t_h
        #  return false
        #elseif curr_l>curr_t_l
        #  nl = max(curr_l, curr_t_l)
        #  nh = min(curr_h, curr_t_h)
        #  msize = nh-nl
        #  osize = curr_h-curr_l
        #  ratio = msize/osize
        #  return rand()>ratio
        #else
        #  return false
        #end


end




#################### Auxilliary Functions ####################

function nutsepsilon(x::T, logfgrad::Function, delta::Float64)  where T<:Node
  n = size(x)[1] - 1
  dir = bitrand(n)
  e1 = nutsepsilon_sub(x, logfgrad, delta, dir)
  #dir = .!dir
  #e2 = nutsepsilon_sub(x, logfgrad, delta, dir)
  e = e1
  #e = 0.0015
  println("e ", e)

  return e
end

function nutsepsilon_sub(x::Node, logfgrad::Function, delta::Float64, dir::BitArray)

  x0 = deepcopy(x)
  n = size(x)[1] - 1

  _, r0, logf0, grad0,_ = refraction(x0, randn(n), 1.0, zeros(n), 0.0, logfgrad, delta, dir, n)

  x0 = deepcopy(x)
  epsilon = 1.0
  _, rprime, logfprime, gradprime,_ = refraction(x0, r0, 1.0 ,grad0,  epsilon, logfgrad, delta, dir,n)

  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  pm = 2 * (prob > 0.5) - 1
  while prob^pm > 0.5^pm
    epsilon *= 2.0^pm
    x0 = deepcopy(x)
    _, rprime, logfprime, _ ,_ = refraction(x0, r0, 1.0, grad0, epsilon, logfgrad, delta, dir, n)

    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  end
  #x = x0
  println("eps ",epsilon)
  epsilon
end
