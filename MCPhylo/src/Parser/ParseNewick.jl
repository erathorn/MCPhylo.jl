

"""
    ParseNewick(filename::String)

This function parses a Newick file
"""
function ParseNewick()
    # for simplicity's sake I split everything in tiny methods
end


# TODO: function to load Newick file and get the string
# Question: should we go with filename ot file-path? just some julia thing?
function load_newick(filename::String)
    open(filename, "r") do file
        global content = readlines(file)
    end
    content[1]
end

# TODO:check if the brackets match etc
function is_valid_newick_string
end

# TODO: create a dataframe, similar to the nexus file
function createNewickdf
end

println(load_newick((pwd()*"\\newick_test_one.txt")))
