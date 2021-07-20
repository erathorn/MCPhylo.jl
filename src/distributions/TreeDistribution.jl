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


"""
    UniformLength()

Fallback struct that is used when no length distribution is given 

"""
struct UniformLength end


"""
    LengthDistribution(length_distr::Union{CompoundDirichlet, exponentialBL, UniformLength}

Wrapper struct for the length distribution. Contains a CompoundDirichlet, exponentialBL or
a uniform length.
"""
mutable struct LengthDistribution
    length_distr::Union{CompoundDirichlet, exponentialBL, UniformLength}
end


"""
    TopologyDistribution(
        constraint_dict::Dict{Symbol, 
                              Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}
                              } where S <: AbstractString)

Wrapper struct for the topology distribution. Contains a dictionary with constraints.
Empty constructor results in an empty dictionary as a fallback.
"""


mutable struct TopologyDistribution
    constraint_dict::Dict{Symbol, 
                          Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}
                         } where S <: AbstractString

    TopologyDistribution() = 
        new(Dict{Symbol, Union{Vector{Vector{AbstractString}}, Vector{Tuple{Vector{AbstractString}, Vector{AbstractString}}}}}()) 
end
    

"""
    TreeDistribution(length_distr::LengthDistribution, topology_distr::TopologyDistribution)

Wrapper struct for the TreeDistribution. Contains a length distribution and a topology
distribution. Can be constructed with no arguments, one of the distributions, or both of 
them. Fallbacks are the UniformLength for the length distribution & and an empty constraint
dictionary for the topology distribution.
"""
mutable struct TreeDistribution <: ContinuousUnivariateDistribution
    length_distr::LengthDistribution
    topology_distr::TopologyDistribution

    TreeDistribution(l::Union{CompoundDirichlet, exponentialBL}, c::Dict) =
        new(LengthDistribution(l), TopologyDistribution(c))
    TreeDistribution(l::Union{CompoundDirichlet, exponentialBL}) =
        new(LengthDistribution(l), TopologyDistribution())
    TreeDistribution(c::Dict) =
        new(LengthDistribution(UniformLength()), TopologyDistribution(c))
    TreeDistribution() = 
        new(LengthDistribution(UniformLength()), TopologyDistribution())
end


function insupport(d::TreeDistribution, x::T) where T <: GeneralNode
    insupport(d.length_distr, x) && insupport(d.topology_distr, x)
end

function logpdf(d::TreeDistribution, x::T) where T <: GeneralNode
    logpdf(d.length_distr.length_distr, x) + logpdf(d.topology_distr, x)
end

function gradlogpdf(d::TreeDistribution, x::T) where T <: GeneralNode
    rl1, rl2 = gradlogpdf(d.length_distr.length_distr, x)
    rt1, rt2 = gradlogpdf(d.topology_distr, x)
    return rl1+rt1, rl2+rt2
end

# placeholder, can be made more elegant using the topology_distr
function relistlength(d::TreeDistribution, x::AbstractArray)
    n = length(x)
    (Array(x), n)
  end