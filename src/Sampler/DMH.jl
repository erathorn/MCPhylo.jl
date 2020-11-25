#################### Double Metropolis-Hastings Sampler ####################

#################### Types and Constructors ####################

mutable struct DMHTune <: SamplerTune
  logf::Union{Function, Missing}
  pseudolog::Union{Function, Missing}
  condlike::Union{Function, Missing}
  datakeys::Vector{Symbol}
  m::Int64
  scale::Float64

  DMHTune() = new()

  function DMHTune(f::Union{Function, Missing}, ps, cl, m::Int64, s::Float64)
    new(f, ps, cl, Vector{Symbol}(), m, s)
  end
  function DMHTune(f::Union{Function, Missing}, ps, cl, dk::Vector{Symbol}, m::Int64, s::Float64)
    new(f, ps, cl, dk, m, s)
  end
end

const DMHVariate = SamplerVariate{DMHTune}

DMHTune(x::Vector, f, ps, cl, dk, m, s) = DMHTune(f, ps, cl, dk, m, s)

#################### Sampler Constructor ####################

function DMH(params::ElementOrVector{Symbol}, m::Int64, scale::Float64)
	println("This")
  samplerfx = function(model::Model, block::Integer)
	println("here")
	tune = gettune(model, block)
	params = asvec(params)
	println(params)
	targets = keys(model, :target, params)
	tune.datakeys = targets
    block = SamplingBlock(model, block, true)
	println(tune.datakeys)

    f = x -> logpdf!(block, x)
	fp = (x, y) -> pseudologpdf!(block, x, y)
	cl = (x, args...) -> conditional_likelihood!(block, x, args...)
    v = SamplerVariate(block, f, fp, cl, targets, m, scale)

    sample!(v::DMHVariate, model)

    relist(block, v)
  end
  Sampler(params, samplerfx, DMHTune())
end

sample!(v::DMHVariate, model) = DMH_sample!(v, v.tune, model)

function DMH_sample!(v::DMHVariate, tune, model)
	# 1. propose theta prime
	# 2. Generate the auxillary variable using theta prime

	obsdata = unlist(model[tune.datakeys[1]])
	obsdata = reshape(obsdata, size(getindex(model, tune.datakeys[1])))
	nfeatures, nlang = size(obsdata)
	@assert v.tune.m >= nlang

	# this way of generating theta_prime from the current values of theta
	# takes care of the transition probability from theta_prime to theta and vice versa
	# the values equal and will cancel out.
	θ_prime = v + v.tune.scale .* rand(Normal(0.0, 1.0), length(v))

	# calculate logpdf values
	lf_xt = tune.logf(v.value)
	lf_xtp = tune.logf(θ_prime)

	# store old value for possible future reference
	θ = v.value

	# prematurely assign, to allow computations to go through properly
	v[:] = θ_prime

	# Sample new pseudo observations based on the θ_prime
	y = inner_sampler(v, deepcopy(obsdata))

	# calculate logpdfs based on pseudo observations

	lf_ytp = tune.pseudolog(v.value, y)
	lf_yt = tune.pseudolog(θ, y)

	# calculate acceptance probability
	r = exp((lf_yt + lf_xtp) - (lf_xt + lf_ytp))

	# RWM acceptance logic is a bit reversed here. Above the new value is
	# prematurely assigned to allow computations in the inner sampler to go through.
	# If the sample is accepted nothing needs to be done anymore, otherwise the
	# old value will be reassigend.
	if rand() > r
		# sample is rejected, so use the old value.
		v[:] = θ
	end
	v
end

function inner_sampler(v::DMHVariate, X)


	nfeatures, nlang = size(X)
	counter = 0
	while true
		random_language_idx = shuffle(1:nlang)
		random_feature_idx = shuffle(1:nfeatures)
		for i in random_language_idx
			for f in random_feature_idx
				feature_vals = sort!(unique(skipmissing(X[f,:])))
				if ismissing(X[i])
					#missing_probs = cond_sum.(Ref(samples), i, f, Ref(spatial_sums),
					#Ref(ling_sums), Ref(v_params), Ref(h_params), Ref(u_params))
					probs = v.tune.condlike(v, i, f)
					missing_probs ./= sum(missing_probs)
					new_x = rand(Categorical(probs))#sample(feature_vals, StatsBase.weights(missing_probs))
					X[f, i] = new_x
				else
					#probs = cond_sum.(Ref(samples), i, f, Ref(spatial_sums), Ref(ling_sums),
					#	Ref(v_params), Ref(h_params), Ref(u_params))
					probs = v.tune.condlike(v, i, f)
					probs ./= sum(probs)
					new_x = rand(Categorical(probs))#sample(feature_vals, StatsBase.weights(probs))
					X[f, i] = new_x
				end
			end
			counter += 1
			if counter > v.tune.m
				return X
			end
		end
	end
end


function cond_sum(X::Array{Union{Missing,Int64},2}, l::Int64, f::Int64,
	spatial_sums::Array{Float64,3}, ling_sums::Array{Float64,3}, v_params, h_params, u_params)
	nvals = maximum(skipmissing(X[f,:]))
	prob_kth_val = Vector{Float64}(undef, nvals)
	for k in 1:nvals
		p = v_params[l] * spatial_sums[f,l,k] + h_params[l] * ling_sums[f,l,k] + u_params[l]
		prob_kth_val[k] = exp(p)
	end
	return prob_kth_val
end
