

mutable struct treevariate{N} <: ArrayVariate{N}
    value::Array{Float64, N}
    symbol::Symbol
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::Mamba.DistributionStruct

end

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
    value = Array{Float64}(undef, fill(0, d)...)
    fx, src = Mamba.modelfxsrc(Mamba.depfxargs, f)
    s = TreeVariate(value, :nothing, Int[], fx, src, Symbol[],
                      Mamba.NullUnivariateDistribution())
    setmonitor!(s, monitor[1], monitor[2])
end

function setmonitor!(d::TreeVariate, monitor::Bool, target)
    println("here")

    d.file = target
    d
end
