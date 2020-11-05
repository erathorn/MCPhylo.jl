#################### Double Metropolis-Hastings Sampler ####################

#################### Types and Constructors ####################

mutable struct DMHTune <: SamplerTune
  outer<:SamplerVariate
  inner<:SamplerVariate

  DMHTune() = new()
  function DMHTune(outer::S, inner::T) where {S<:SamplerVariate, T<:SamplerVariate}
    typeof(outer) == DMHVariate && throw(ArgumentError($outer " cannot be of type " $typeof(outer)))
    typeof(inner) == DMHVariate && throw(ArgumentError($inner " cannot be of type " $typeof(inner)))
    new(outer, inner)
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
    v = SamplerVariate(block, f, NullFunction(); args...)

    DMH_sample!(v::DMHVariate, f)

    relist(block, v)
  end
  Sampler(params, samplerfx, DMHTune())
end

function DMH_sample!(v::DMHVariate, lf::Function, pslf::Function)
  # 1. propose theta prime
  # 2. Generate the auxillary variable using theta prime

  # this way of generating theta_prime from the current values of theta
  # takes care of the transition probability from theta_prime to theta and vice versa
  # the values equal and will cancel out.
  θ_prime = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))

	lf_xt = lf(v.value)
	lf_xtp = lf(θ_prime)

	### get the original data into the sampler
	y = inner_sampler(v)

	lf_yt = pslf(v.value, y)
	lf_ytp = pslf(θ_prime, y)

	r = exp((lf_yt + lf_xtp) - (lf_xt + lf_ytp))

	if rand() < r
		v[:] = x
	end
	v
end

function m_inner_sampler(v::DMHVariate)

	samples = similar()
end



function inner_sampler(x, v, h, u, spatial_graph, ling_graph, m)
	nlang, n = size(spatial_graph)
	feature_vals = sort!(unique(x))
	samples = Array{Int64,1}()
	idx_list = Vector(1:length(x))
	while m != 0 # is this the best way to limit the number of iterations?
		for i in x
			i = rand(idx_list) # primitive randomisation - problem: returns a vector
			if ismissing(x[i]) # no missing values yet, and this part needs improvement;
				# To do: add the neighbour and global majority methods
				missing_probs = Array{Float64,1}()
				missing_probs_norm = Array{Float64,1}()
				for k in feature_vals
					p_unnorm = cond_prob(x, k, i, spatial_graph, ling_graph, dummy_v, dummy_h, dummy_u)
					push!(missing_probs, p_unnorm)
				end
				for p in missing_probs
					p_norm = p/sum(missing_probs)
					push!(missing_probs_norm, p_norm)
				end
				new_x = sample(feature_vals, StatsBase.weights(missing_probs_norm))
				push!(samples, new_x) # with random iteration, how do I make sure samples[i] corresponds to the right x[i]?
			else
				probs = Array{Float64,1}()
				probs_norm = Array{Float64,1}() # inefficient?
				for k in feature_vals
					p_unnorm = cond_prob(x, k, i, spatial_graph, ling_graph, dummy_v, dummy_h, dummy_u)
					push!(probs, p_unnorm)
				end
				for p in probs
					p_norm = p/sum(probs)
					push!(probs_norm, p_norm)
				end
				new_x = sample(feature_vals, StatsBase.weights(probs_norm))
				push!(samples, new_x)
				filter!(x->x != i, idx_list)
			end
			m -= 1
		end
		return samples
	end
end
