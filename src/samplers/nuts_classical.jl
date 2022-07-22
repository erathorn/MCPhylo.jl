


function nuts_sub!(v::Sampler{NUTSTune{classic_nuts, F}, T}, epsilon::R, logfgrad::Function) where {T<: AbstractArray{<: Real}, F, R<:Real}
  n = length(v)
  x = deepcopy(v.value)
  logf, grad = logfgrad(x)
  r = randn(n)
  
  logp0 = logf - 0.5 * turbo_dot(r, r)
  logu0 = logp0 + log(rand())
  xminus = xplus = deepcopy(x)
  rminus = rplus = deepcopy(r)
  gradminus = gradplus = deepcopy(grad)
  j = 0
  n = 1
  s = true
  while s && j < v.tune.tree_depth
    pm = 2 * (rand() > 0.5) - 1
    if pm == -1
      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0)
    else
      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha =
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0)
    end
    if sprime && rand() < nprime / n
      v[:] = xprime
    end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus)
    v.tune.stepsizeadapter.metro_acc_prob = alpha / nalpha
  end
  v
end


function leapfrog(x::Vector{Float64}, r::Vector{Float64}, grad::Vector{Float64},
                  epsilon::Real, logfgrad::Function)
  r += (0.5 * epsilon) * grad
  x += epsilon * r
  logf, grad = logfgrad(x)
  r += (0.5 * epsilon) * grad
  x, r, logf, grad
end


function buildtree(x::Vector{Float64}, r::Vector{Float64},
                   grad::Vector{Float64}, pm::Integer, j::Integer,
                   epsilon::R, logfgrad::Function, logp0::R, logu0::R) where R<:Real
  if j == 0
    xprime, rprime, logfprime, gradprime = leapfrog(x, r, grad, pm * epsilon,
                                                    logfgrad)
    logpprime = logfprime - 0.5 * turbo_dot(rprime, rprime)
    nprime = Int(logu0 < logpprime)
    sprime = logu0 < logpprime + 1000.0
    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    alphaprime = isnan(alphaprime) ? 0.0 : alphaprime
    nalphaprime = 1
  else
    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0)
    if sprime
      if pm == -1
        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
      end
      if rand() < nprime2 / (nprime + nprime2)
        xprime = xprime2
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
    end
  end
  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime
end


function nouturn(xminus::Vector{Float64}, xplus::Vector{Float64},
                 rminus::Vector{Float64}, rplus::Vector{Float64})
  xdiff = xplus - xminus
  turbo_dot(xdiff, rminus) >= 0 && turbo_dot(xdiff, rplus) >= 0
end
