

function make_tree_with_data_cu(filename::String, dialect::AbstractString="nexus",
                             gap::Union{Missing, AbstractString}=missing,
                             miss::Union{Missing,AbstractString}=missing,
                             header::Bool=false)
	new_tree, my_df = make_tree_with_data(filename, dialect, gap, miss, header)
 	my_df = CuArray{Float64}(my_df)

    return new_tree, my_df
end # function make_tree_with_data


function make_tree_with_data(filename::String, dialect::AbstractString="nexus",
                             gap::Union{Missing, AbstractString}=missing,
                             miss::Union{Missing,AbstractString}=missing,
                             header::Bool=false;
							 binary::Bool=true)
    # get all the information from the input file
	if lowercase(dialect) == "nexus"
		n_tax, n_sites, gap, miss, symbols, df = ParseNexus(filename)
	elseif lowercase(dialect) == "csv"
		ismissing(gap) || throw("Please specify the gap symbol for a CSV file")
		ismissing(miss) || throw("Please specify the missing symbol for a CSV file")
		n_tax, n_sites, symbols, df = ParseCSV(filename, header)
	end
	# create random tree
	new_tree = create_tree_from_leaves(df[!,:Language], n_sites)

	n_nodes = length(post_order(new_tree))
	n_states = length(symbols)
	my_df = Array{Float64}(undef, n_states, n_sites, n_nodes)


	# iterate through the data frame and get the node information
	for row in eachrow(df)
		#data_vec = zeros(Float64, (2, n_sites))
		mn = find_by_name(new_tree, row.Language)
		mind = mn.num
		for (ind, i) in enumerate(row.Data)
			ent = string(i)
			index = findfirst(x->x==ent, symbols)
			if index == nothing
				if ent == gap || ent == miss
					my_df[:, ind, mind] .= 1.0
				else
					println(typeof(ent))
					println(typeof.(symbols))
					throw("unknown symbol $ent, $symbols")
				end
			else
				my_df[:,ind,mind] .= 0.0
				my_df[index,ind,mind] = 1.0
			end # if
		end # for
		#node = find_by_name(new_tree, row.Language)
		#node.data = log.(data_vec)
	end # for


    return new_tree, my_df
end # function make_tree_with_data
