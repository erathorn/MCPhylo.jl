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

# const extensiontests = [
#   "newunivardist",
#   "newmultivardist"
# ]



println("Running tests:")
#
@testset "All tests" begin


for t in samplertests
    @runtest "samplers/" t
  end

@testset "mcmctests" begin
   for t in mcmctests
     @runtest "mcmc/" t
   end
end

end
#all sets test set
