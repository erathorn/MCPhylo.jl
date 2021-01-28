

function parse_and_number(treestring::S)::Node where S <: AbstractString
    p_tree2 = parsing_newick_string(string(treestring))
    set_binary!(p_tree2)
    number_nodes!(p_tree2)
    p_tree2
end


function nodnums2names(bps::Vector{Tuple}, tree::Node)::Vector{Tuple{Set{String}, Set{String}}}
    treebps = Vector{Tuple{Set{String}, Set{String}}}()
    for bipart in bps
        bipartition = []
        for part in bipart
            myset = Set{String}()
            for num in part
                no = find_num(tree, num)
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
function ASDSF(filename1::S, filename2::S, freq::Int64=1)::Vector{Float64} where S <: AbstractString
    splitsQueue = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueue_1 = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    splitsQueue_2 = Accumulator{Tuple{Set{String}, Set{String}}, Int64}()
    ASDF_vals = zeros(Int(countlines(filename1)/freq))
    all_keys = Set{Tuple{Set{String},Set{String}}}()
    run = 1
    p = Progress(length(ASDF_vals))
    for ((ind1, tree_string_1), (ind2, tree_string_2)) in zip(enumerate(eachline(filename1)), enumerate(eachline(filename2)))
        if mod(ind1,freq) == 0
            tree1 = parse_and_number(tree_string_1)
            tree2 = parse_and_number(tree_string_2)
            bp_t1 = get_bipartitions(tree1)
            bp_t2 = get_bipartitions(tree2)
            named_bp_1 = nodnums2names(bp_t1, tree1)
            named_bp_2 = nodnums2names(bp_t2, tree2)
            cm1_d = countmap(named_bp_1)
            cm2_d = countmap(named_bp_2)
            outer = 0.0
            new_splits = union(keys(cm1_d), keys(cm2_d))
            all_keys = union(all_keys, new_splits)
            for split in new_splits
                inc!(splitsQueue, split)
            end
            for split in all_keys
                inner = 0.0
                if haskey(cm1_d, split)
                    inner += 0.5
                    inc!(splitsQueue_1, split)
                end
                if haskey(cm2_d, split)
                    inner += 0.5
                    inc!(splitsQueue_2, split)
                end
                if splitsQueue_1[split]/run > 0.1 || splitsQueue_2[split]/run > 0.1
                    fre = splitsQueue[split]/run
                    inner += (inner - fre)^2
                    outer += sqrt(inner)
                end
            end
            ASDF_vals[run] = outer / length(splitsQueue)
            run += 1
            next!(p)
        end
    end
    ASDF_vals
end
