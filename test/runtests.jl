using Distributed
using MCPhylo
using Random, Test
include("utils.jl")

const distributionstests = [
    "phylodist"
]

const likelihoodtests = [
    "rates",
    "substitutionmodels",
    "felsenstein"
]

const mcmctests = [
    "discretediag",
    "modelbuilding"
#   "readcoda"
]

const modeltests = [
    "test_dependent"
]

const outputtests = [
    "test_chains"
]

const parsertests = [
    "csv",
    "nexus",
    "parser"
]

const samplertests = [
    #"abc",
    "linear_model_uvp",
    "linear_model_mvp"
]

println("Running tests:")


@testset "All tests" begin

    @testset "distributions tests" begin
        for t in distributionstests
            @runtest "distributions/" t
        end
    end
    
    @testset "likelihood tests" begin
        for t in likelihoodtests
            @runtest "likelihood/" t
        end
    end

    @testset "mcmc tests" begin
        for t in mcmctests
            @runtest "mcmc/" t
        end
    end

    @testset "model tests" begin
        for t in modeltests
            @runtest "model/" t
        end
    end

    @testset "output tests" begin
        for t in outputtests
            @runtest "output/" t
        end
    end

    @testset "parser tests" begin
        for t in parsertests
            @runtest "parser/" t
        end
    end

    @testset "sampler tests" begin
        for t in samplertests 
            @runtest "samplers/" t
        end
    end
end
