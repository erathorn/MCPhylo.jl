#################### LengthDistributions ####################

"""
    exponentialBL(scale::Float64) <: ContinuousUnivariateDistribution

This structure implememts an exponential prior on the branch lengths of a tree.
"""
struct exponentialBL <: ContinuousUnivariateDistribution
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
struct CompoundDirichlet <: ContinuousUnivariateDistribution
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
    UniformBranchLength()

Fallback struct that is used when no length distribution is given 

"""
struct UniformBranchLength <: ContinuousUnivariateDistribution end

const LengthDistribution = Union{CompoundDirichlet, exponentialBL, UniformBranchLength}

#################### TopologyDistributions ####################

"""
    UniformConstrained(
        constraint_dict::Dict{Symbol, 
                              Union{Vector{Vector{S}}, Vector{Tuple{Vector{S}, Vector{S}}}}
                              } where S <: AbstractString)

Wrapper struct for the topology distribution. Contains a dictionary with constraints.
Empty constructor results in an empty dictionary as a fallback.
"""
struct UniformConstrained <: ContinuousUnivariateDistribution
    constraint_dict::ConstraintDict

    UniformConstrained() = 
        new(ConstraintDict()) 
    UniformConstrained(x::Dict) = 
        new(ConstraintDict(x)) 
end

"""
    UniformTopology()

Fallback struct that is used when no topology distribution is given.

"""
struct UniformTopology <: ContinuousUnivariateDistribution end

const TopologyDistribution = Union{UniformConstrained, UniformTopology}


"""
    TreeDistribution(length_distr::LengthDistribution, topology_distr::TopologyDistribution)

Wrapper struct for the TreeDistribution. Contains a length distribution and a topology
distribution. Can be constructed with no arguments, one of the distributions, or both of 
them. Fallbacks are UniformBranchLength for the length distribution & and UniformTopology
for the topology distribution.
"""
struct TreeDistribution <: ContinuousUnivariateDistribution
    length_distr::LengthDistribution
    topology_distr::TopologyDistribution


    TreeDistribution(l::LengthDistribution) =
        new(l, UniformTopology())
    TreeDistribution(l::LengthDistribution, c::Dict) =
        new(l, UniformConstrained(ConstraintDict(c)))
    TreeDistribution(c::ConstraintDict) =
        new(UniformBranchLength(), UniformConstrained(c))
    TreeDistribution(c::Dict) =
        new(UniformBranchLength(), UniformConstrained(ConstraintDict(c)))
    TreeDistribution() = 
        new(UniformBranchLength(), UniformTopology())
end


function insupport(d::TreeDistribution, x::GeneralNode)
    insupport(d.length_distr, x) && insupport(d.topology_distr, x)
end

function logpdf(d::TreeDistribution, x::GeneralNode, transform::Bool)
    logpdf(d.length_distr, x) + logpdf(d.topology_distr, x)
end

function gradlogpdf(d::TreeDistribution, x::GeneralNode)
    rl1, rl2 = gradlogpdf(d.length_distr, x)
    rt1, rt2 = gradlogpdf(d.topology_distr, x)
    return rl1+rt1, rl2 .+ rt2
end

# placeholder, can be made more elegant using the topology_distr
function relistlength(d::TreeDistribution, x::AbstractArray)
    x[1], 1
  end