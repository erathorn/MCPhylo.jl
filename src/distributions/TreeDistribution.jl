

"""
 TODO:
    - DistributionType for length and topology distribution
    - proper constructors

"""
mutable struct TreeDistribution <: ContinuousUnivariateDistribution
    length_distr
    topology_distr
end

function insupport(d::TreeDistribution, x::T) where T <: GeneralNode
    insupport(d.length_distr, x) && insupport(d.topology_distr, x)
end

function logpdf(d::TreeDistribution, x::T) where T <: GeneralNode
    logpdf(d.length_distr, x) + logpdf(d.topology_distr, x)
end

function gradlogpdf(d::TreeDistribution, x::T) where T <: GeneralNode
    rl1, rl2 = gradlogpdf(d.length_distr, x)
    rt1, rt2 = gradlogpdf(d.topology_distr, x)
    return rl1+rt1, rl2+rt2
end

# placeholder, can be made more elegant using the topology_distr
function relistlength(d::TreeDistribution, x::AbstractArray)
    n = length(x)
    (Array(x), n)
  end