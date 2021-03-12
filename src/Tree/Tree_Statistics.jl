function parse_and_number(treestring::S)::Node where S <: AbstractString
    p_tree2 = parsing_newick_string(string(treestring))
    set_binary!(p_tree2)
    number_nodes!(p_tree2)
    p_tree2
end # parse_and_number


function nodnums2names(bps::Vector{Tuple}, tree::Node
                      )::Vector{Tuple{Set{String}, Set{String}}}

    treebps = Vector{Tuple{Set{String}, Set{String}}}()
    for bipart in bps
        bipartition = []
        for part in bipart
            myset = Set{String}()
            for num in part
                no = find_num(tree, num)
                if no.nchild == 0
                    push!(myset, no.name)
                end # if
            end # for
            push!(bipartition, myset)
        end # for
        push!(treebps, Tuple(bipartition))
    end # for
    treebps
end # nodnums2names


"""
    ASDSF(args...; freq::Int64=1, check_leafsets::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two files
containing newick representations of trees.
"""
function ASDSF(args...; freq::Int64=1, check_leafsets::Bool=true
              )::Vector{Float64}

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
            check_leafsets && check_leafsets(trees)
            bps = [get_bipartitions(tree) for tree in trees]
            named_bps = [nodnums2names(bp, tree) for (bp, tree) in zip(bps, trees)]
            cmds = [countmap(named_bp) for named_bp in named_bps]
            outer = 0.0
            new_splits = union([keys(cmds[t]) for t in 1:length(trees)]...)
            all_keys = union(all_keys, new_splits)
            for split in new_splits
                inc!(splitsQueue, split)
            end # for
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


"""
    ASDSF(args...; freq::Int64=1, check_leafsets::Bool=true)::Vector{Float64}

Calculate the average standard deviation of split frequencies for two files
containing newick representations of trees.
"""
function ASDSF(args...; freq::Int64=1, check_leafsets::Bool=true
              )::Vector{Float64}

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
            check_leafsets && check_leafsets(trees)
            bps = [get_bipartitions(tree) for tree in trees]
            named_bps = [nodnums2names(bp, tree) for (bp, tree) in zip(bps, trees)]
            cmds = [countmap(named_bp) for named_bp in named_bps]
            outer = 0.0
            new_splits = union([keys(cmds[t]) for t in 1:length(trees)]...)
            all_keys = union(all_keys, new_splits)
            for split in new_splits
                inc!(splitsQueue, split)
            end # for
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
