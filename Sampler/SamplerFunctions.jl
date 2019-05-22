module SamplerFunctions

using Markdown
using Distributions

export SampleVector!, SampleVectorAdjusted!

"""
        SampleVector(parameters::Vector{Float64} ind::Integer, sd::Float64, lower::Float64, upper::Float64)

This function samples a new value for the value at position `ind` in the
parameter vector. The vector is changed inplace and the logarithm of the
Metropolis-Hastings ratio is returned. The parameters `lower` and `upper`
are the sampling boundaries.
The new values are sampled using a normal distribution which is centered at the
old value and the standard deviation of this distribution is specified in `sd`.
"""
function SampleVector!(parameters::Vector{Float64}, ind::Integer, sd::Float64, lower::Float64, upper::Float64)
        my_value::Float64 = parameters[ind]
        sampling_dist = Distributions.Normal(my_value, sd)
        while true
                new_value = rand(sampling_dist)
                if lower<= new_value <= upper
                        parameters[ind] = new_value
                        break
                end # if
        end # while

        return Distributions.logpdf(Distributions.Normal(parameters[ind],sd))-Distributions.logpdf(Distributions.Normal(my_value,sd))
end # function SampleVector

"""
        SampleVectorAdjusted!(parameters::Vector{Float64}, ind::Integer, sd::Float64, lower::Float64, upper::Float64, ov_sum::Float64)

This function samples a new value for the value at position `ind` in the
parameter vector. The vector is changed inplace and the logarithm of the
Metropolis-Hastings ratio is returned. The parameters `lower` and `upper`
are the sampling boundaries.
The new values are sampled using a normal distribution which is centered at the
old value and the standard deviation of this distribution is specified in `sd`.
Additionally, this function makes sure that the sum of the values in `parameters`
remains unchanged. The sum of the parameter values is specified in `ov_sum`.
"""
function SampleVectorAdjusted!(parameters::Vector{Float64}, ind::Integer, sd::Float64, lower::Float64, upper::Float64, ov_sum::Float64)
        my_value::Float64 = parameters[ind]
        sampling_dist = Distributions.Normal(my_value, sd)
        sample::Bool = true
        while sample
                new_value = rand(sampling_dist)
                if lower<= new_value <= upper
                        scale = (ov_sum-new_value)/(ov_sum-new_value)
                        parameters[ind] = new_value
                        sample = false
                        for i in 1:len(parameters)
                                if i != ind
                                        parameters[i] *= scale
                                end # if
                        end # for
                        break
                end # if
        end # while
        return Distributions.logpdf(Distributions.Normal(new,sd))-Distributions.logpdf(Distributions.Normal(my_value,sd))
end # function

end  # modul SamplerFunctions
