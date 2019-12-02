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

  function PNUTSTune(x::Vector{Node}, epsilon::Float64, logfgrad::Union{Function, Missing};
                    target::Real=0.6)
    new(logfgrad, false, 0.0, epsilon, 1.0, 0.05, 0.0, 0.75, 0, NaN, 0, 10.0,0.003,
        target,0)
  end
end

PNUTSTune(x::Vector{Node}, logfgrad::Function, ::NullFunction, delta::Float64=0.003; args...) =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)

PNUTSTune(x::Vector{Node}, logfgrad::Function, delta::Float64; args...) =
  PNUTSTune(x, nutsepsilon(x[1], logfgrad, delta), logfgrad; args...)


  NUTSTune(x::Vector{Node}, epsilon::Float64; args...) =
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

function mlogpdfgrad!(block::SamplingBlock, x::Node, sz::Int64, ll::Bool=false, gr::Bool=false)
  grad = Vector{Float64}(undef, sz)
  lp = 0.0
  if ll
    lp = logpdf!(block, x)



  end
  if gr
    grad = gradf!(block, x, :provided)



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
  lor = 0.5 > rand()
  x = deepcopy(mt)
  currU, _ = logfgrad(x, nl, true, false)
  r = randn(nl)
  org_r = deepcopy(r)
  currH = currU-0.5*dot(r)
  x, r , logf, grad, nni = refraction(mt, r, zeros(nl), epsilon, logfgrad,delta, lor, nl)

  probH = logf-0.5*dot(r)
  ratio = currH - probH
  if  ratio >= min(0, log(rand()))
    logp0 = logf - 0.5 * dot(r)
    logu0 = logp0 + log(rand())
    xminus = xplus = x
    rminus = rplus = r
    gradminus = gradplus = grad
  else
    logp0 = currU - 0.5 * dot(r)
    logu0 = logp0 + log(rand())
    xminus = xplus = mt
    rminus = rplus = org_r
    gradminus = gradplus = zeros(nl)
  end

  j = 0
  n = 1
  s = true

  while s


    pm = 2 * (rand() > 0.5) - 1
    if pm == -1
      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha, nni1, lpp = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0, delta, false, nl)
    else
      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha, nni1, lpp =
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0, delta, true, nl)
    end
    if sprime && rand() < nprime / n
      #if min(1, exp(lpp-logp0)) > rand()
        v.value[1]= xprime
    #end
  end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus, grad, epsilon, logfgrad, delta, nl, j)
    v.tune.alpha, v.tune.nalpha = alpha, nalpha
    nni += nni1
  end

  v
end

function scale_fac(x::T, delta::T) where T<:Real
  x < delta ? x/delta : 1.0
end


function refraction(v::Node, r::Vector{Float64},
                    grad::Vector{Float64}, epsilon::Float64, logfgrad::Function,
                    delta::Float64, lor::Bool, sz::Int64)

    ref_grad = -grad
    blenvec = zeros(sz)
    get_branchlength_vector(v, blenvec)
    fac = scale_fac.(blenvec, delta)


    r = @. r - epsilon * 0.5 * (ref_grad * fac)

    tmpB = @. blenvec + (epsilon * r)

    nni = 0

    if minimum(tmpB) <= 0
        v, tmpB, r, nni = ref_NNI(v, tmpB, r, epsilon, blenvec, delta, logfgrad, lor, sz)

    end

    set_branchlength_vector!(v, molify(tmpB, delta))
    logf, grad = logfgrad(v, sz, true, true)
    ref_grad = -grad
    get_branchlength_vector(v, blenvec)
    fac = scale_fac.(blenvec, delta)
    r = @. r - epsilon * 0.5 * (ref_grad * fac)

    return v, r, logf, grad, nni
end


function ref_NNI(v::Node, tmpB::Vector{Float64}, r::Vector{Float64}, epsilon::Float64, rv::Vector{Float64},
                 delta::Float64, logfgrad::Function, lor::Bool, sz::Int64)

  intext = internal_external(v)
  t = 0.0
  nni = 0
  ## r = probM
  while minimum(tmpB)<=0.0

     timelist = tmpB./abs.(r)
     ref_index = argmin(timelist)

     temp = epsilon-t+timelist[ref_index]
     rv = @. rv + temp * r

     r[ref_index] *= -1.0

     if intext[ref_index] == 1
       blv =molify(rv, delta)
       set_branchlength_vector!(v, blv)
       U_before_nni, _ = logfgrad(v, sz, true, false) # still with molified branch length
       U_before_nni *= -1
       v_copy = deepcopy(v)
       tmp_NNI_made = NNI!(v_copy, ref_index, lor)
       if tmp_NNI_made != 0

            U_after_nni, _ = logfgrad(v_copy, sz, true, false)
            U_after_nni *= -1
            delta_U::Float64 = 2.0*(U_after_nni - U_before_nni)
            my_v::Float64 = r[ref_index]^2
            if my_v >= delta_U

              r[ref_index] = sqrt(my_v - delta_U)

              v = v_copy
              nni += tmp_NNI_made
            end # if my_v
        end #if NNI
      end #non leave
      t = epsilon+timelist[ref_index]

      tmpB = @. rv + (epsilon-t) * r

  end #while

  v, tmpB, r, nni
end



function buildtree(x::Node, r::Vector{Float64},
                   grad::Vector{Float64}, pm::Integer, j::Integer,
                   epsilon::Float64, logfgrad::Function, logp0::Real, logu0::Real,
                   delta::Float64, lor::Bool, sz::Int64)


  if j == 0


    xprime, rprime, logfprime, gradprime, nni = refraction(x, r, pm.*grad, epsilon,
                                          logfgrad, delta, lor, sz)

    logpprime = logfprime - 0.5 * dot(rprime)
    nprime = Int(logu0 < logpprime && min(1, exp(logpprime - logp0)) > rand())
    sprime = logu0 < logpprime + 1000.0
    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    nalphaprime = 1
  else

    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime, nni, logpprime = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0, delta, lor, sz)
    if sprime
      if pm == -1

        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 , nni, logpprime= buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta, lor, sz)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2, nni, logpprime = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0, delta, lor, sz)
      end
      if rand() < nprime2 / (nprime + nprime2)
        if min(1, exp(logpprime - logp0)) > rand()
          xprime = xprime2
        end
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus, grad, epsilon, logfgrad, delta, sz, j)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
    end
  end

  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime, nni, logpprime
end


function nouturn(xminus::Node, xplus::Node,
                rminus::Vector{Float64}, rplus::Vector{Float64}, grad::Vector{Float64},
                epsilon::Float64, logfgrad::Function, delta::Float64, sz::Int64, j::Int64)

        if j > 7
          return false
        end
        rf0 = RF(xminus, xplus)
        xminus_bar,_,_,_,_ = refraction(deepcopy(xminus), rminus, -grad, epsilon, logfgrad, delta, false, sz)
        xplus_bar,_,_,_,_ = refraction(deepcopy(xplus), rplus, grad, epsilon, logfgrad, delta, true, sz)
        rf1 = RF(xminus_bar, xplus_bar)
        if rf1 == rf0
          xminusone = zeros(sz)
          xplusone = zeros(sz)
          xminus_barone = zeros(sz)
          xplus_barone = zeros(sz)
          get_branchlength_vector(xminus_bar, xminus_barone)
          get_branchlength_vector(xplus_bar, xplus_barone)
          get_branchlength_vector(xminus, xminusone)
          get_branchlength_vector(xplus, xplusone)

          xdiff1 = xplusone-xminusone
          xdiff2 = xplus_barone-xminus_barone

          return dot(xdiff1) >= dot(xdiff2)

        else

          return rf1 > rf0
        end

end




#################### Auxilliary Functions ####################

function nutsepsilon(x::Node, logfgrad::Function, delta::Float64)
  e1 = nutsepsilon_sub(x, logfgrad, delta, true)
  e2 = nutsepsilon_sub(x, logfgrad, delta, false)
  e = (e1+e2)/2.0
  println("e ", e)

  return e
end

function nutsepsilon_sub(x::Node, logfgrad::Function, delta::Float64, dir::Bool)

  x0 = deepcopy(x)
  n = size(x)[1] - 1

  _, r0, logf0, grad0,_ = refraction(x, randn(n), zeros(n), 0.0, logfgrad, delta, dir, n)

  x = deepcopy(x0)
  epsilon = 1.0
  _, rprime, logfprime, gradprime,_ = refraction(x, r0, grad0, epsilon, logfgrad, delta, dir,n)

  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  pm = 2 * (prob > 0.5) - 1
  while prob^pm > 0.5^pm
    epsilon *= 2.0^pm
    x = deepcopy(x0)
    _, rprime, logfprime, _ ,_ = refraction(x, r0, grad0, epsilon, logfgrad, delta, dir, n)

    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))

  end
  x = x0
  println("eps ",epsilon)
  epsilon
end
