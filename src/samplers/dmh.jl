#################### Double Metropolis-Hastings Sampler ####################

#################### Types and Constructors ####################

mutable struct DMHTune <: SamplerTune
  outer<:SamplerVariate
  inner<:SamplerVariate
  seed_data<:AbstractArray
  m<:Int64

  DMHTune() = new()
  function DMHTune(outer::S, inner::T, seed_data::A, m::Int64) where {S<:SamplerVariate, T<:SamplerVariate}
    typeof(outer) == DMHVariate && throw(ArgumentError($outer " cannot be of type " $typeof(outer)))
    typeof(inner) == DMHVariate && throw(ArgumentError($inner " cannot be of type " $typeof(inner)))
	@assert m >= length(seed_data)
    new(outer, inner, seed_data, m)
  end
end

const DMHVariate = SamplerVariate{DMHTune}


#################### Sampler Constructor ####################

function DMH(params::ElementOrVector{Symbol}, outer::Symbol,
             inner::Symbol)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
	targets = keys(model, :target, block)
    f = x -> logpdf!(block, x)
	fp = (x, y) -> pseudologpdf!(block, x, y)
	cl = x -> conditional_likelihood(block, x)
    v = SamplerVariate(block, f, NullFunction(); args...)

    DMH_sample!(v::DMHVariate, f, fp, cl)

    relist(block, v)
  end
  Sampler(params, samplerfx, DMHTune())
end

function DMH_sample!(v::DMHVariate, lf::Function, pslf::Function, cl::Function)
	# 1. propose theta prime
	# 2. Generate the auxillary variable using theta prime

	# this way of generating theta_prime from the current values of theta
	# takes care of the transition probability from theta_prime to theta and vice versa
	# the values equal and will cancel out.
	θ_prime = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))

	# calculate logpdf values
	lf_xt = lf(v.value)
	lf_xtp = lf(θ_prime)

	# store old value for possible future reference
	θ = v.value

	# prematurely assign, to allow computations to go through properly
	v[:] = θ_prime

	# Sample new pseudo observations based on the θ_prime
	y = inner_sampler(v, cl)

	# calculate logpdfs based on pseudo observations
	lf_ytp = pslf(v.value, y)
	lf_yt = pslf(θ, y)

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

function inner_sampler(v::DMHTune, cond_prob::Function)
	#X::Array{Int64,1}, v_params, h_params,
	#u_params, spatial_sums::Array{Float64,2}, ling_sums::Array{Float64,2}, m::Int64)
	m = v.m
	X = v.seed_data

	feature_vals = sort!(unique(X))
	N = length(X)
	samples = deepcopy(X)
	idx_list = Vector(1:N)
	counter = 0
	while true
		random_index_order = shuffle(1:N)
		for i in random_index_order
			if ismissing(X[i]) # no missing values yet, and this part needs improvement;
				# To do: add the neighbour and global majority methods
				throw("We should not end up here")
				missing_probs = cond_prob.(Ref(samples), i, feature_vals, Ref(spatial_sums),
					Ref(ling_sums), Ref(v_params), Ref(h_params), Ref(u_params))
				missing_probs ./= sum(missing_probs)
				new_x = sample(feature_vals, StatsBase.weights(missing_probs))
				samples[i] = new_x
			end
				probs = cond_prob.(Ref(samples), i, feature_vals, Ref(spatial_sums), Ref(ling_sums),
					Ref(v_params), Ref(h_params), Ref(u_params))
				probs ./= sum(probs)
				new_x = sample(feature_vals, StatsBase.weights(probs))
				samples[i] = new_x
			end
			counter += 1
			if counter > m
				return samples
			end
		end
	end
end
