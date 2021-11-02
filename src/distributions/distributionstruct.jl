#################### DistributionStruct ####################

#################### Base Methods ####################

dims(d::DistributionStruct) = size(d)

function dims(D::Array{MultivariateDistribution})
  size(D)..., mapreduce(length, max, D)
end


#################### List Fallbacks ####################

unlist(d::Distribution, x::AbstractArray) = vec(x)

unlist_sub(d::Distribution, x::AbstractArray) = unlist(d, x)

unlist_sub(d::UnivariateDistribution, X::AbstractArray) = vec(X)

unlist_sub(D::Array{UnivariateDistribution}, X::AbstractArray) = vec(X)

function unlist_sub(D::Array{MultivariateDistribution}, X::AbstractArray)
  y = similar(X, length(X))
  offset = 0
  for sub in CartesianIndices(size(D))
    n = length(D[sub])
    inds = 1:n
    y[offset .+ inds] = X[sub, inds]
    offset += n
  end
  resize!(y, offset)
end

function relistlength(d::UnivariateDistribution, x::AbstractArray)
  (x[1], 1)
end

function relistlength(d::MultivariateDistribution, x::AbstractArray)
  n = length(d)
  value = x[1:n]
  (value, n)
end

function relistlength(d::MatrixDistribution, x::AbstractArray)
  n = length(d)
  value = reshape(x[1:n], size(d))
  (value, n)
end

relistlength_sub(d::Distribution, s::AbstractStochastic, x::AbstractArray) = relistlength(d, x)

function relistlength_sub(d::Distribution, s::Stochastic{T}, x::T) where T<:GeneralNode
  relistlength(d, x)
end

function relistlength_sub(d::Distribution, s::Stochastic{T}, x::AbstractArray{T}) where T<:GeneralNode
  relistlength(d, x[1])
end

relistlength(d::UnivariateDistribution, x::T) where T<:GeneralNode = (x, 1)

function relistlength_sub(d::Union{Array{UnivariateDistribution}, UnivariateDistribution}, 
                          s::Stochastic{<:AbstractArray{<:Real, N} where N}, X::AbstractArray)
  n = length(s)
  value = reshape(X[1:n], size(s))
  (value, n)
end

function relistlength_sub(D::Array{MultivariateDistribution}, s::Stochastic{<:AbstractArray}, X::AbstractArray)
  Y = similar(X, size(s))
  offset = 0
  for sub in CartesianIndices(size(D))
    n = length(D[sub])
    inds = 1:n
    Y[sub, inds] = X[offset .+ inds]
    offset += n
  end
  (Y, offset)
end


#################### Link Fallbacks ####################

# link(d::Distribution, x) = x

link_sub(d::Distribution, x) = link(d, x)

function link_sub(d::PDMatDistribution, x)
  res = mylink(d, x)

  res
end


# function link_sub(d::UnivariateDistribution, X::AbstractArray)
#   Y = similar(X, Float64)
#   map!(x -> link_sub(d, x), Y, X)
# end

# function link_sub(D::Array{UnivariateDistribution}, X::AbstractArray)
#   Y = similar(X, Float64)
#   map!(i -> link_sub(D[i], X[i]), Y, 1:length(D))
# end

# function link_sub(D::Array{MultivariateDistribution}, X::AbstractArray)
#   Y = similar(X, Float64)
#   for sub in CartesianIndices(size(D))
#     d = D[sub]
#     inds = 1:length(d)
#     Y[sub, inds] = link_sub(d, X[sub, inds])
#   end
#   Y
# end


# invlink(d::Distribution, x) = x

invlink_sub(d::Distribution, x) = invlink(d, x)

function invlink_sub(d::PDMatDistribution, x) 
  U, _ = relistlength(d, unlist(d, myinvlink(d, x)))
  U
end



# function invlink_sub(d::UnivariateDistribution, X::AbstractArray)
#   #Y = similar(X, Float64)
#   map(x -> invlink_sub(d, x), X)
# end

# function invlink_sub(D::Array{UnivariateDistribution}, X::AbstractArray)
#   Y = similar(X, Float64)
#   map!(i -> invlink_sub(D[i], X[i]), Y, 1:length(D))
# end

# function invlink_sub(D::Array{MultivariateDistribution}, X::AbstractArray)
#   Y = similar(X, Float64)
#   for sub in CartesianIndices(size(D))
#     d = D[sub]
#     inds = 1:length(d)
#     Y[sub, inds] = invlink_sub(d, X[sub, inds])
#   end
#   Y
# end


logcond(d, x, args...) = logcond(d, x, args...)

#################### Logpdf Fallbacks ####################

logpdf(d::Distribution, x, transform::Bool) = logpdf(d, x)

function logpdf_sub(d::Distribution, x, transform::Bool)
  insupport(d, x) ? logpdf(d, x, transform) : -Inf
end

function logpdf_sub(d::DiscreteMatrixDistribution, X::AbstractArray, transform::Bool)
  logpdf(d, X)
end

function logpdf_sub(d::UnivariateDistribution, X::AbstractArray,
                    transform::Bool)
  lp = sum([logpdf_sub(d, X[i], transform) for i in 1:length(X)])
  lp
end

function logpdf_sub(D::Array{UnivariateDistribution}, X::AbstractArray,
                    transform::Bool)
  @inbounds lp = sum([logpdf_sub(D[i], X[i], transform) for i in 1:length(D)])
  lp
end

function logpdf_sub(D::Array{MultivariateDistribution}, X::AbstractArray,
                    transform::Bool)
  
  @inbounds lp = sum([logpdf_sub(D[i], vec(X[i, :]), transform) for i in 1:size(D,1)])
  lp
end

function gradlogpdf_sub(D::Array{UnivariateDistribution}, X::AbstractArray)
    if length(D) > 1
      throw(DimensionMismatch("To many Distributions"))
    end
    grad = gradlogpdf_sub(D[1], X[1])
    ifelse.(isfinite.(grad), grad, 0.0)
end

function gradlogpdf_sub(d::Distribution, x)
  gradlogpdf(d, x)
end

#################### Rand Fallbacks ####################

rand_sub(d::Distribution, x) = rand(d)

rand_sub(d::UnivariateDistribution, X::AbstractArray) = rand(d, size(X))

rand_sub(D::Array{UnivariateDistribution}, X::AbstractArray) = map(rand, D)

function rand_sub(D::Array{MultivariateDistribution}, X::AbstractArray)
  Y = fill(NaN, size(X))
  for sub in CartesianIndices(size(D))
    d = D[sub]
    Y[sub, 1:length(d)] = rand(d)
  end
  Y
end