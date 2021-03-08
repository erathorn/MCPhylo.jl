

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

"""
this is the function I've been using, as it actually returns results I can turn into a confusion matrix;
the other function just replaces unknown values according to the result of FelsensteinFunction()
"""
function fill_in_the_blanks_test(filename::String,ratesfile::String, tree::String)

	ntax,nchar,_,_,df = ParseNexus(filename)

	ratefile_contents = []
	open(ratesfile, "r") do file
	  for x in readlines(file)
	    push!(ratefile_contents,split(x,"\t"))
	  end #for
	end #open
	for x in ratefile_contents
		if x[1] == "pi(0)"
			global pi_ = parse(Float64,x[6])
		end #if
		if x[1] == "alpha"
			global rates = parse(Float64,x[6])
		end #if
	end #for


	langnames = df[1]
	numbers = df[2]

	if tree == ""
		tree, _= make_tree_with_data(filename, replace_missing=false)
	else
		tree = parsing_newick_string(tree)
	end #ifelse

	number_nodes!(tree)

	# standardizes node.data
	for node_to_edit in post_order(tree)
		node_to_edit.data = Array{Float64}(undef,2,nchar)
	end #for

	#assigns data values from dataframe to nodes in tree; easier to work with
	for x in eachindex(numbers)
		try
			global cur_node = find_by_name(tree,langnames[x])
		catch
			println("Dataframe/tree mismatch; make sure data and tree are properly formatted!")
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


	global results = []
    global pairwise_results = []
	global already_tested = []

	for trial in range(1,stop=10000)
		#selects random language/index in node.data that isn't unknown OR already used
		while true
			global randlang = rand(collect(range(1,ntax)))
			global randnum = rand(collect(range(1,nchar)))
			if numbers[randlang][randnum] == '?' || (randlang,randnum) in already_tested
				continue
			else
				break
			end #ifelse
		end #while
		push!(already_tested, (randlang,randnum))
		node_of_interest = find_by_name(tree, langnames[randlang])
		global trueval = node_of_interest.data[1,randnum]
		#temporarily assign 3.0/unknown values to the index in node.data we're interested in testing
		node_of_interest.data[1,randnum] = 3.0
		node_of_interest.data[2,randnum] = 3.0
		testval = assign_value(node_of_interest.name,tree,randnum,ntax,nchar,pi_,rates)
		if trueval == testval
			push!(results,true)
		else
			push!(results,false)
		end #ifelse
		push!(pairwise_results, [trueval,testval])
		#reassign the true value of node.data
		node_of_interest.data[1,randnum] = trueval
		node_of_interest.data[2,randnum] = Float64(!Bool(trueval))
	end #for
	return results,pairwise_results, tree
end #fill_in_the_blanks_test



function fill_in_the_blanks(filename::String, ratesfile::String, tree::String)

	ntax,nchar,_,_,df = ParseNexus(filename)

	langnames = df[1]
	numbers = df[2]

	ratefile_contents = []
	open(ratesfile, "r") do file
	  for x in readlines(file)
	    push!(ratefile_contents,split(x,"\t"))
	  end #for
	end #open
	for x in ratefile_contents
		if x[1] == "pi(0)"
			global pi_ = parse(Float64,x[6])
		end #if
		if x[1] == "alpha"
			global rates = parse(Float64,x[6])
		end #if
	end #for

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
			println("Dataframe/tree mismatch; make sure data and tree are properly formatted!")
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
				retval = assign_value(x.name,tree,y,ntax,nchar,pi_,rates)
				x.data[1,y] = Float64(retval)
				x.data[2,y] = Float64(!retval)
			end #if
		end #for
	end #for
	return tree
end #function fill_in_the_blanks
"""
This function actually calls FelsensteinFunction() and is used in both fill_in_the_blank() functions
"""
function assign_value(langname::String,tree::Node,IoI::Int64,ntax::Int64,nchar::Int64,pi_::Float64,rates::Float64)
	treecopy = deepcopy(tree)
	#global already_removed = []
	#this for loop handles pruning; I know it's not perfectly implemented, but it works and I didn't want to change it until I figured out the
	#issue with false negatives
	# for node in post_order(treecopy)
	# 	#"if node has an unknown value in the same index as the node we're inferring a value for, isn't the root,
	# 	#isn't the same node we're inferring a value for, is a leaf, and hasn't already been removed..."
	# 	if node.data[1,IoI]  == 3.0 && !node.root && node.name != langname && isempty(node.children) && !(node.num in already_removed)
	# 		while true
	# 			#this while loop prunes superfluous internal nodes
	# 			mom = node.mother
	# 			push!(already_removed,node.num)
	# 			remove_child!(node.mother,node)
	# 			if isempty(mom.children)
	# 				node = mom
	# 			else
	# 				break
	# 			end #ifelse
	# 		end #while
	# 	end #if
	# end #for

	for node in get_leaves(treecopy)
		if node.name != langname && node.data[1,IoI] == 3.0
			mom = node.mother
			# println( "FOLLOWING NODE PURGED: ")
			# println(node.name)
			# println("OFFENDING DATAPOINT: ", node.data[1,IoI], "MATCHES NoI ", find_by_name(tree,langname).name, "'s ", find_by_name(tree,langname).data[1,IoI], "AT POSITION ", IoI)
			remove_child!(mom,node)
			while length(mom.children) == 1 && !(mom.root)
				grandma = mom.mother
				child = mom.children[1]
				remove_child!(mom,child)
				add_child!(grandma,child)
				remove_child!(grandma, mom)
				mom = grandma
			end #while
			while isempty(mom.children)
				oldmom = mom
				mom = mom.mother
				# println("BLV VALUE:")
				# println(oldmom.blv)
				remove_child!(mom,oldmom)
			end #while
		end #if
	end #for
	number_nodes!(treecopy)
	global pruned_tree = deepcopy(treecopy)
	treecopy = reroot(treecopy,langname)
	global rerooted_tree = deepcopy(treecopy)
	global tree_collection = [newick(tree),newick(pruned_tree),newick(rerooted_tree)]

	number_nodes!(treecopy)
	blv = get_branchlength_vector(treecopy)
	dataplaceholder = Array{Float64}(undef,1,1,1)
	final_rates = discrete_gamma_rates(rates,rates,4)
	felsenstein_results = []
	for ind_rate in final_rates
		# println(tree_collection[1], "\n\n", tree_collection[2], "\n\n", tree_collection[3])
		FelsensteinFunction(post_order(treecopy), pi_, ind_rate, dataplaceholder, nchar, blv)

		# if treecopy.data[1,IoI] == NaN
		# 	println(treecopy.name)
		# 	println(IoI)
		# 	println(treecopy.data[1,:])
		# 	println(treecopy.data[2,:])
		# end #if
		push!(felsenstein_results,treecopy.data[1,IoI])
		treecopy.data[1,IoI] = 3.0
	end #for
	if sum(felsenstein_results) >= 0.5
		return true
	else
		return false
	end #ifelse
	return false
end #assign_value
