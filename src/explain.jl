include("./MCPhylo.jl")
using .MCPhylo
using DataStructures
using ProgressMeter
using StatsBase

asdsf = ASDSF_old("./src/tree1.nwk", "./src/tree3.nwk")

trees = MCPhylo.ParseNewick("./doc/Tree/Drav_mytrees_1.nwk")

"""
plot1 = Plots.plot(trees[1])
plot2 = Plots.plot(trees[1], treetype=:fan, msc=:blue, mc=:yellow, lc=:white,
           bg=:black, tipfont=(7, :lightgreen))
"""

data = rand(Normal(0,1), 5000)

my_data=Dict(:data=>data)

model = Model(
    data = Stochastic(1, (μ, σ) -> Normal(μ, σ), false),
       μ = Stochastic(()->Normal(),true),
       σ = Stochastic(()->Exponential(1), true)
)

inits = [Dict(:data => data,
            :μ => randn(),
            :σ => rand()),
       Dict(:data => data,
           :μ => randn(),
           :σ => rand())]

samplers = [NUTS(:μ),
           Slice(:σ, 0.1)]

setsamplers!(model, samplers)

sim = mcmc(model, my_data, inits, 1000, burnin=500, thin=5, chains=2, trees=true)

# default "inner" layout puts plots in a row
pv = plot(sim, [:mean])
# "inner" layout can be manipulated, but usually size has to be adjusted as well
pv = plot(sim, [:mean], layout=(3, 1), size=(800,1500))
# throws an error, as it should for contour (when only one variable is selected)
pv = plot(sim, [:contour], vars=["likelihood"])
# gives a warning for contourplot but shows the other ptypes
pv = plot(sim, [:contour, :density, :mean], vars=["likelihood"], fuse=true)
# specific plot variables are passed successfully
pv = plot(sim, [:autocor, :contour, :density, :mean, :trace],
           maxlag=10, bins=20, trim=(0.1, 0.9), legend=true)
# demonstrate the customizable "outer" layout
pv = plot(sim, [:autocor, :bar, :contour, :mixeddensity, :mean, :trace],
           fuse=true, fLayout=(2,2), fsize=(2750, 2500), linecolor=:match)
# barplot works
pv = plot(sim, [:bar], linecolor=:match, legend=:true, filename="blub.pdf")
# use savefig to save as file; no draw function needed
savefig("test.pdf")





function parse_and_number(treestring::S)::Node where S <: AbstractString
    p_tree2 = MCPhylo.parsing_newick_string(string(treestring))
    MCPhylo.set_binary!(p_tree2)
    MCPhylo.number_nodes!(p_tree2)
    p_tree2
end


function nodnums2names(bps::Vector{Tuple}, tree::Node)::Vector{Tuple{Set{String}, Set{String}}}
    treebps = Vector{Tuple{Set{String}, Set{String}}}()
    for bipart in bps
        bipartition = []
        for part in bipart
            myset = Set{String}()
            for num in part
                no = MCPhylo.find_num(tree, num)
                if no.nchild == 0
                    push!(myset, no.name)
                end
            end
            push!(bipartition, myset)
        end
        push!(treebps, Tuple(bipartition))
    end
    treebps
end


"""
    ASDSF(filename1::S, filename2::S, freq::Int64=1)::Vector{Float64} where S <: AbstractString

Calculate the average standard deviation of split frequencies for two files containing
newick representations of trees.
"""
function ASDSF(args...; freq::Int64=1)::Vector{Float64} where S <: AbstractString
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueues = Vector{Accumulator{Tuple{Set{String}, Set{String}}, Int64}}()
    for arg in args
        push!(splitsQueues, Accumulator{Tuple{Set{String}, Set{String}}, Int64}())
    end
    ASDF_vals = zeros(Int(countlines(args[1])/freq))
    all_keys = Set{Tuple{Set{String},Set{String}}}()
    run = 1
    p = Progress(length(ASDF_vals))
    for (i, lines) in enumerate(zip([eachline(arg) for arg in args]...))
        if mod(i, freq) == 0
            trees = [parse_and_number(tree_string) for tree_string in lines]
            check_leafsets(trees)
            bps = [MCPhylo.get_bipartitions(tree) for tree in trees]
            named_bps = [nodnums2names(bp, tree) for (bp, tree) in zip(bps, trees)]
            cmds = [countmap(named_bp) for named_bp in named_bps]
            outer = 0.0
            new_splits = union([keys(cmds[t]) for t in 1:length(trees)]...)
            all_keys = union(all_keys, new_splits)
            for split in new_splits
                inc!(splitsQueue, split)
            end
            for split in all_keys
                inner = 0.0
                for (t, tree) in enumerate(trees)
                    if haskey(cmds[t], split)
                        inner += 1 / length(trees)
                        inc!(splitsQueues[t], split)
                    end # if
                end # for
                if any([x[split] / run > 0.1 for x in splitsQueues])
                    fre = splitsQueue[split] / run
                    inner += (inner - fre) ^ 2
                    outer += sqrt(inner)
                end # if
            end # for
            ASDF_vals[run] = outer / length(splitsQueue)
            run += 1
            ProgressMeter.next!(p)
        end # if
    end # for
    ASDF_vals
end # ASDSF
