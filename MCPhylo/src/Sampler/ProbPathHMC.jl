



mutable struct ProbPathHMCTune <: SamplerTune
    n_leap::Int64
    stepsz::Float64

    ProbPathHMCTune() = new()

    ProbPathHMCTune(n_leap::Int64, stpesz::Float64) = new(n_leap::Int64, stpesz::Float64)
end # mutable struct

const ProbPathVariate = SamplerVariate{ProbPathHMCTune}


function ProbPathHMCSampler(params, pargs...; dtype::Symbol=:forward)

    samplerfx = function(model::Model, block::Integer)

        #println(model[model.samplers[block].params[1]])

        #block = SamplingBlock(model, block, true)
        v = ProbPathHMCTune(pargs...)
        println(size(model.nodes[:data]))
        #t = ProbPathVariate(model[model.samplers[block].params[1]], v)
        #SamplerVariate(v)
        #t = ProbPathVariate(v)
        sample!(v, model.nodes[:mtree], model.nodes[:data], model.nodes[:mypi], model.nodes[:mtree].distr)
        #println(v)
        #relist(block, v)
    end # function samplerfx
    Sampler(params, samplerfx, ProbPathHMCTune())
end

function sample!(v::ProbPathHMCTune, tree ,data, mypi, distr)
    println("here")
    n_c = size(data)[2]
    my_sample!(tree, data, v.n_leap, v.stepsz, mypi, n_c, distr)
end


function Stochastic(d::Integer, f::Function, monitor::Tuple)
    value = Array{Float64}(undef, fill(0, 2)...)
    fx, src = modelfxsrc(depfxargs, f)
    s = TreeVariate(value, :nothing, Int[], fx, src, Symbol[],
                      NullUnivariateDistribution(), "ab")
    setmonitor!(s, monitor[1], monitor[2])
end

function setmonitor!(d::TreeVariate, monitor::Bool, target::AbstractString)
    println("here")

    d.file = target
    d.monitor = [1,2]
    d
end

"""
    setinits!(d::TreeVariate, m::model, x::Array)

documentation
"""
function setinits!(d::TreeVariate, m::Model, x::Array)
    d.value = x
    d.distr = eval(m)
    setmonitor!(d, d.monitor)
end # function
