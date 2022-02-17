using Distributed
using MCPhylo
using Random, Test
include("utils.jl")


const samplertests = [
  #"abc",
  "linear_model_uvp",
  "linear_model_mvp"
]

const mcmctests = [
   "discretediag",
   "modelbuilding"
#   "readcoda"
]


const parsertests = [
  "csv",
  "nexus"
]


# const extensiontests = [
#   "newunivardist",
#   "newmultivardist"
# ]

println("Running tests:")


@testset "All tests" begin

  @testset "sampler tests" begin
    for t in samplertests 
      @runtest "samplers/" t
    end
  end

  @testset "mcmc tests" begin
    for t in mcmctests
      @runtest "mcmc/" t
    end
  end

  @testset "parser tests" begin
    for t in parsertests
      @runtest "parser/" t
    end
  end

end
