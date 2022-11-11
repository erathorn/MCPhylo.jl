

"""
	make_tree_with_data(filename::String, dialect::AbstractString="nexus",
							 gap::Union{Missing, AbstractString}=missing,
							 miss::Union{Missing,AbstractString}=missing,
							 header::Bool=false; binary::Bool=true)::Tuple{GeneralNode, Array{Float64}}

General parsing function; user specifies type of file to parse.
Returns root Node of tree derived from data, as well as a DataFrame.

* `filename` : Name of file to parse.
* `dialect` : Specifies filetype to parse; currently accepts EITHER "nexus" or "csv" as inputs.
* `gap` : specifies gap symbol for a CSV file; if "dialect" == "csv", this input is required, can be ommitted otherwise.
* `miss` : specifies missing symbol for a CSV file; if "dialect" == "csv", this input is required, can be ommitted otherwise.
* `header` : Boolean required for parsing CSV file, specifies if file contains a header; input not required when parsing a NEXUS file.
* `binary` : Boolean, specifies whether a binary or nonbinary tree should be created for the given file. Defaults to true/binary tree.
"""
function make_tree_with_data(
    filename::String,
    dialect::AbstractString = "nexus",
    gap::Union{Missing,AbstractString} = missing,
    miss::Union{Missing,AbstractString} = missing,
    header::Bool = false;
    binary::Bool = true,
)::Tuple{GeneralNode,Array{Float64}}

    # get all the information from the input file
    if lowercase(dialect) == "nexus"
        n_tax, n_sites, gap, miss, symbols, df, langs = ParseNexus(filename)
    elseif lowercase(dialect) == "csv"
        ismissing(gap) &&
            throw(ArgumentError("Please specify the gap symbol for a CSV file"))
        ismissing(miss) &&
            throw(ArgumentError("Please specify the missing symbol for a CSV file"))
        n_tax, n_sites, gap, miss, symbols, df, langs =
            ParseCSV(filename, gap, miss, header)
    end

    new_tree = create_tree_from_leaves(langs, binary)

    my_df = datafortree(df, langs, new_tree, symbols, gap, miss)

    return new_tree, my_df
end # function make_tree_with_data


function datafortree(
    df::Matrix{Char},
    leave_names::Vector{A},
    tree::T,
    symbols::Vector{C},
    gap::A,
    miss::A;
    log_space::Bool = false,
)::Array{
    Float64,
} where {T<:GeneralNode,A<:AbstractString,B<:AbstractString,C<:AbstractString}
    n_nodes = length(post_order(tree))
    n_states = length(symbols)
    n_sites = size(df, 2)
    my_df = Array{Float64}(undef, n_states, n_sites, n_nodes)


    # iterate through the data frame and get the node information
    for (ind, row) in enumerate(eachrow(df))
        mn = find_by_name(tree, leave_names[ind])
        mind = mn.num
        for (ind, i) in enumerate(row)
            ent = string(i)
            index = findfirst(x -> x == ent, symbols)
            if index === nothing
                if ent == gap || ent == miss
                    my_df[:, ind, mind] .= log_space ? log(1.0) : 1.0
                else
                    throw("unknown symbol $ent, $symbols")
                end
            else
                my_df[:, ind, mind] .= log_space ? log(0.0) : 0.0
                my_df[index, ind, mind] = log_space ? log(1.0) : 1.0
            end # if
        end # for
    end # for
    my_df
end
