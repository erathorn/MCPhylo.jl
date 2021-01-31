

function make_tree_with_data_cu(filename::String, dialect::AbstractString="nexus",
                             gap::Union{Missing, AbstractString}=missing,
                             miss::Union{Missing,AbstractString}=missing,
                             header::Bool=false;
							 replace_missing::Bool=true
							 )
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
			 elseif i == '?'
				 if replace_missing == true
					 my_df[1, ind, mind] = 1.0
					 my_df[2, ind, mind] = 1.0
				 else
					 my_df[1,ind,mind] = 3.0
					 my_df[2,ind,mind] = 3.0
				 end #if
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
							 binary::Bool=true,
							 replace_missing::Bool=true
							 )
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
			elseif i == '?'
				if replace_missing == true
					my_df[1, ind, mind] = 1.0
					my_df[2, ind, mind] = 1.0
				else
					my_df[1,ind,mind] = 3.0
					my_df[2,ind,mind] = 3.0
				end #if

            else
                my_df[1, ind, mind] = 1.0
                my_df[2, ind, mind] = 1.0

            end # if
        end # for
    end # for

    return new_tree, my_df
end # function make_tree_with_data

function fill_in_the_blanks(filename::String, tree::String="")

	ntax,nchar,_,_,df = ParseNexus(filename)

	langnames = df[1]
	numbers = df[2]

	if tree == ""
		tree, _= make_tree_with_data(filename, replace_missing=false)
	else
		tree = parsing_newick_string(tree)
	end #ifelse
	number_nodes!(tree)
	for node_to_edit in post_order(tree)
		node_to_edit.data = Array{Float64}(undef,2,nchar)
	end #for

	for x in eachindex(numbers)
		try
			global cur_node = find_by_name(tree,langnames[x])
		catch
			continue
		end #try
		for y in eachindex(numbers[1])
			if numbers[x][y] != '?'
				cur_node.data[1,y] = Float64(Bool(parse(Int,numbers[x][y])))
				cur_node.data[2,y] = Float64(!Bool(parse(Int,numbers[x][y])))
			else
				cur_node.data[1,y] = 3.0
				cur_node.data[2,y] = 3.0
			end #ifelse
		end #for
	end #for


	for x in post_order(tree)
		for y in eachindex(x.data[1,:])
			if x.data[1,y] == 3.0
				retval = assign_value(x.name,tree,y,ntax,nchar)
				x.data[1,y] = Float64(retval)
				x.data[2,y] = Float64(!retval)
				# println(x.data[1,y])
				# println("\n\n\n")
			end #if
		end #for
	end #for
	return tree,retdf
end #function fill_in_the_blanks

function assign_value(langname::String,tree::Node,IoI::Int64,ntax::Int64,nchar::Int64)
	global uncopied_leaves = get_leaves(tree)
	for x in uncopied_leaves
		if x.name == "no_name"
			throw("MEST UP")
		end
	end
	treecopy = deepcopy(tree)
	global copied_leaves = get_leaves(treecopy)
	for x in copied_leaves
		if x.name == "no_name"
			throw("MEST UP")
		end
	end

	#global nodes2prune = Node[]
	treecopy = reroot(treecopy,langname)
	global rerooted_leaves = get_leaves(treecopy)
	for x in rerooted_leaves
		if x.name == "no_name"
			throw("MEST UP")
		end
	end

	for node in post_order(treecopy)
		if node.data[1,IoI]  == 3.0 && !node.root && node.name != langname && isempty(node.children)
			while true
				mom = node.mother
				remove_child!(node.mother,node)
				if isempty(mom.children)
					node = mom
				else
					break
				end #ifelse
			end #while
			#push!(nodes2prune,node)
		end #if
	end #for
	global pruned_tree = get_leaves(treecopy)
	for x in pruned_tree
		if x.name == "no_name"
			throw("MEST UP")
		end
	end

	# if !isempty(nodes2prune)
	# 	for removeme in nodes2prune
	# 		remove_child!(removeme.mother,removeme)
	# 	end #for
	# end #if
	number_nodes!(treecopy)
	blv = get_branchlength_vector(treecopy)
	dataplaceholder = Array{Float64}(undef,1,1,1)
	thing = treecopy.data[1,:]
	 lengthodata = treecopy.data[1,:]


	 pi_= (count(i->(i==1.0),thing))/nchar
	try
		FelsensteinFunction(post_order(treecopy), pi_, [1.0], dataplaceholder, nchar, blv)
	catch
		println("log error")
	end
	# println(treecopy.data[1,IoI])
	if treecopy.data[1,IoI] >= 0.5
		return true
	else
		return false
	end #ifelse
	return false





	# rerooted_tree = reroot(tree,lang_of_interest
	# )
	# blv = get_branchlength_vector(rerooted_tree)
	# println("do i get here?")
	# println("rerooted_tree: ", "\t", typeof(post_order(rerooted_tree)))
	# println("pi_: ", "\t", typeof(1))
	# println("rates: ", "\t", typeof([1.0]))
	# data = [1.0,1.0,1.0]
	# println("data: ", "\t", typeof(data))
	# println("n_char: ", "\t", typeof(n_char))
	# println("blv: ", "\t", typeof(blv))
	# likelihood = FelsensteinFunction(post_order(rerooted_tree), 1, [1.0], data, n_char, blv)
	# if likelihood > 0.5
	# 	dataframe[dim1,dim2,dim3] = 1.0
	# else
	# 	dataframe[dim1,dim2,dim3] = 0.0
	# end # if
end #assign_value
