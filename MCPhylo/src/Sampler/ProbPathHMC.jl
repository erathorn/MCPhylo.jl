



mutable struct ProbPathHMCTune <: SamplerTune
    n_leap::Float64
    stepsz::Float64

    ProbPathHMCTune() = new()

    ProbPathHMCTune(n_leap::Float64, stpesz::Float64) = new(n_leap::Float64, stpesz::Float64)
end # mutable struct

const ProbPathVariate = SamplerVariate{ProbPathHMCTune}


function ProbPathHMCSampler(params, pargs...; dtype::Symbol=:forward)

    samplerfx = function(model::Model, block::Integer)

        #println(model[model.samplers[block].params[1]])

        #block = SamplingBlock(model, block, true)
        v = ProbPathHMCTune(pargs...)

        #t = ProbPathVariate(model[model.samplers[block].params[1]], v)
        #SamplerVariate(v)
        #t = ProbPathVariate(v)
        sample!(v, model)
        #println(v)
        #relist(block, v)
    end # function samplerfx
    Sampler(params, samplerfx, ProbPathHMCTune())
end

function sample!(v::ProbPathHMCTune, model::Model)
    tree = model.nodes[:mtree]
    mypi = model.nodes[:mypi]
    distr = model.nodes[:mtree].distr
    #n_c = tree.value.data[2]
    nt = my_sample!(tree.value, v.n_leap, v.stepsz, mypi, 3132, distr)
    tree.value = nt
    v
end

macro promote_treevariate(V)
  quote
    Base.promote_rule(::Type{$(esc(V))}, ::Type{T}) where T<:Real = Float64
  end
end
@promote_treevariate TreeVariate

function Stochastic(d::AbstractString, f::Function, monitor::Union{Bool, Vector{Int}}=true)
    value = Node()
    fx, src = modelfxsrc(depfxargs, f)
    s = TreeStochastic(value, :nothing, Int[], fx, src, Symbol[],
                      NullUnivariateDistribution())
    setmonitor!(s, monitor)
end

function setmonitor!(d::TreeStochastic, monitor::Bool)
    d.monitor = [1,2]
    d
end

"""
    setinits!(d::TreeVariate, m::model, x::Array)

documentation
"""
function setinits!(d::TreeStochastic, m::Model, x::Array)
    d.value = x
    d.distr = d.eval(m)
    setmonitor!(d, d.monitor)
end # function


#function relistlength(d::TreeStochastic, v::SubArray, w::Bool)
#
#    ms = size(d.distr)
#    rs = reshape(v, ms)
#    (rs, length(d.distr))
#end

function update!(d::TreeStochastic, m::Model)
    d.distr = d.eval(m)
    d
end

function names(d::TreeStochastic, nodekey::Symbol)
    AbstractString["Tree height", "Tree length"]
end

function unlist(d::TreeStochastic)
    tree_height(d.value), tree_length(d.value)
end

function unlist(s::AbstractStochastic, x::Node, transform::Bool=false)
    s.value
end
