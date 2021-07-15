"""
    exponentialBL(scale::Float64) <: ContinuousUnivariateDistribution
This structure implememts an exponential prior on the branch lengths of a tree.
"""
mutable struct exponentialBL <: ContinuousUnivariateDistribution
    scale::Float64
    constraints::Union{Dict, Missing}

    exponentialBL(scale::Float64) =
        new(scale, missing)
    exponentialBL(scale::Float64, c) =
            new(scale, c)
end


"""
    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64, nterm::Float64)
This structure implememts the CompoundDirichlet distribution described
in Zhang, Rannala and Yang 2012. (DOI:10.1093/sysbio/sys030)
"""
mutable struct CompoundDirichlet <: ContinuousUnivariateDistribution
    alpha::Float64
    a::Float64
    beta::Float64
    c::Float64
    constraints::Union{Dict, Missing}

    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64) =
        new(alpha, a, beta, c, missing)
    CompoundDirichlet(alpha::Float64, a::Float64, beta::Float64, c::Float64, constraints::Dict) =
            new(alpha, a, beta, c, constraints)
end # CompoundDirichlet


mutable struct TopologyDistribution
    constraint_dict::Union{Dict{Symbol, Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}}  where S <: AbstractString, Missing}
end


mutable struct LengthDistribution
    length_distr::Union{CompoundDirichlet, exponentialBL, Missing}
end
    

mutable struct TreeDistribution <: ContinuousUnivariateDistribution
    length_distr::LengthDistribution
    topology_distr::TopologyDistribution

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