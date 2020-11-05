#################### Double Metropolis-Hastings Inner Sampler ####################

#################### Types and Constructors ####################

mutable struct DMHInnerTune <: SamplerTune
  #cond_prob::Function
  m::Int64

  DMHTInnerune() = new()

  #function DMHInnerTune(cond_prob::Function, m::Int64)
  #  new(cond_prob, m)
  #end
	function DMHInnerTune(m::Int64)
    new(m)
  end

end

const DMHVInnerariate = SamplerVariate{DMHInnerTune}


#################### Sampler Constructor ####################

function DMHInner(params::ElementOrVector{Symbol}, m::Int64)

  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    f = (x, l, k)  -> logcond(block, x, l, k)
    v = SamplerVariate(block, f, NullFunction())

    sample!(v::DMHInnerVariate, f)

    relist(block, v)
  end
  Sampler(params, samplerfx, DMHInnerTune())
end

function sample!(v::DMHInnerVariate, logcond::Function)

end

function sample_x(x, v, h, u, spatial_graph, ling_graph, m)
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
