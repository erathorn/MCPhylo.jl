include("./Tree/Tree_Module.jl")
include("./Substitution/SubstitutionMat.jl")
include("./Parser/ParseNexus.jl")
using .Tree_Module
using .SubstitutionMat
using .NexusParser
include("./Likelhood/LikelihoodCalculator.jl")
using .LikelihoodCalculator


this_tree = NexusParser.make_tree_with_data("./local/development.nex")

po = Tree_Module.post_order(this_tree)

#function Samplepi(mpi::Float64, sd::Float64, lower::Float64, upper::Float64)
#        new_value::Float64 = 0
#        sampling_dist = Distributions.TruncatedNormal(mpi, sd, lower, upper)
#        while true
#                new_value = rand(sampling_dist)
#                if lower<= new_value <= upper
#                        break
#                end # if
#        end # while

#        return new_value, Distributions.logpdf(Distributions.TruncatedNormal(new_value,sd, lower, upper), mpi)-
#        Distributions.logpdf(Distributions.TruncatedNormal(mpi,sd, lower, upper), new_value)
#end # function SampleVector






nr = 1000
my_pi = 0.3
cr = LikelihoodCalculator.FelsensteinFunction(po, my_pi)
#cr_arr = []
#mypi_arr = []
#push!(cr_arr, cr)
#push!(mypi_arr, my_pi)
#i = 0
#Samplepi(my_pi, 0.1, 0.0, 1.0)
#while i < nr
#        global my_pi
#        global cr
#        global cr_arr
#        global mypi_arr
#        global i
#        nv, mh = Samplepi(my_pi, 0.1, 0.0, 1.0)

#        nr = LikelihoodCalculator.FelsensteinFunction(po, nv)
#
#        u = log(rand(Distributions.Uniform(0.0, 1.0)))
#        a = min(0.0, (nr-cr)+mh)
#        if u <= a
#                my_pi = nv
#                cr = nr
#        end
#        push!(cr_arr, cr)
#        push!(mypi_arr, my_pi)
#        i += 1
#        println(i)

#end # while
