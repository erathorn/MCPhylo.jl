using Distributed
using MCPhylo
using LinearAlgebra
using FiniteDiff
using Random, Test
include("utils.jl")

const distributionstests = [
    "phylodist",
    "treeconstraints",
    "treedists"
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
    "dependent",
    "variate"
]

const outputtests = [
    "asdsf",
    "chains",
    "gelman"
]

const parsertests = [
    "csv",
    "nexus",
    "parser"
]

const samplertests = [
    #"abc",
    "sampler_construction",
    "linear_model_uvp",
    "linear_model_mvp",
    "tree_samplers"
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
