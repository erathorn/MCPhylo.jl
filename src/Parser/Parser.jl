

function make_tree_with_data_cu(filename::String, dialect::AbstractString="nexus",
                             gap::Union{Missing, AbstractString}=missing,
                             miss::Union{Missing,AbstractString}=missing,
                             header::Bool=false)
		 # get all the information from the input file
     if lowercase(dialect) == "nexus"
         n_tax, nc, gap, miss, df = ParseNexus(filename)
     elseif lowercase(dialect) == "csv"
         ismissing(gap) || throw("Please specify the gap symbol for a CSV file")
         ismissing(miss) || throw("Please specify the missing symbol for a CSV file")
         n_tax, nc, df = ParseCSV(filename, header)
     end
     # create random tree
     new_tree = create_tree_from_leaves_cu(df[!,:Language], nc)

     n_nodes = length(post_order(new_tree))
     #my_df = Array{Float64}(undef, 2, nc, n_nodes)
 	 my_df = Array{Float64}(undef, 2, nc, n_nodes)
     #my_df = ForwardDiff.Dual{Float64}.(my_df, 1)
     #my_df = MArray{Tuple{2,nc,n_nodes}, Float64}(undef)

     # iterate through the data frame and get the node information
     for row in eachrow(df)
         #data_vec = zeros(Float64, (2, nc))
         mn = find_by_name(new_tree, row.Language)
         mind = mn.num
         for (ind, i) in enumerate(row.Data)

             if i == '0'
                 my_df[1,ind,mind] = 1.0
                 my_df[2,ind, mind] = 0.0
             elseif i == '1'
                 my_df[1,ind,mind] = 0.0
                 my_df[2,ind, mind] = 1.0
             else
                 my_df[1, ind, mind] = 1.0
                 my_df[2, ind, mind] = 1.0

             end # if
         end # for
         #node = find_by_name(new_tree, row.Language)
         #node.data = log.(data_vec)
     end # for



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
        n_tax, nc, gap, miss, df = ParseNexus(filename)
    elseif lowercase(dialect) == "csv"
        ismissing(gap) || throw("Please specify the gap symbol for a CSV file")
        ismissing(miss) || throw("Please specify the missing symbol for a CSV file")
        n_tax, nc, df = ParseCSV(filename, header)
    end
    # create random tree
	new_tree = Node()
	if binary
    	new_tree = create_tree_from_leaves_bin(df[!,:Language], nc)
	else
		new_tree = create_tree_from_leaves(df[!,:Language], nc)
	end

    n_nodes = length(post_order(new_tree))
	my_df = Array{Float64}(undef, 2, nc, n_nodes)

    # iterate through the data frame and get the node information
    for row in eachrow(df)

        mn = find_by_name(new_tree, row.Language)
        mind = mn.num
        for (ind, i) in enumerate(row.Data)

            if i == '0'
                my_df[1,ind,mind] = 1.0
                my_df[2,ind, mind] = 0.0
            elseif i == '1'
                my_df[1,ind,mind] = 0.0
                my_df[2,ind, mind] = 1.0
            else
                my_df[1, ind, mind] = 1.0
                my_df[2, ind, mind] = 1.0

            end # if
        end # for
    end # for

    return new_tree, my_df
end # function make_tree_with_data
