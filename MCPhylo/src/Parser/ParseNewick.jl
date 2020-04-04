

"""
    ParseNewick(filename::String)

This function parses a Newick file
"""
function ParseNewick(filename::String)
    # for simplicity's sake I split everything in tiny methods
    content = load_newick(filename)
    if !is_valid_newick_string(content)
        throw("$filename is not a Newick file!")
    #the parse thing


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
function is_valid_newick_string(newick::String)
    # Step one: does the stripped string ends with ';'
    if endswith(strip(newick),";")
        bracket_level = 0
        for char in newick
            if char == '('
                bracket_level += 1
            elseif char == ')'
                    bracket_level -= 1
            end #elseif
        end #for
        if bracket_level != 0
            return false
        end #if

    #end # if (endswith one)
    else
        return false
    end # else
    return true
end

# TODO: create a dataframe, similar to the nexus file
function createNewickdf
end

test_bad_brackets = load_newick((pwd()*"\\newick_test_brackets_fucked.txt"))
test_bad_ending = load_newick((pwd()*"\\newick_test_no_column.txt"))
test_all_gut = load_newick((pwd()*"\\newick.nwk"))
test_new_lines = load_newick((pwd()*"\\newick_new_lines.txt"))
println(test_new_lines)
println(is_valid_newick_string(test_new_lines))
# println(is_valid_newick_string(test_bad_brackets))
# println(is_valid_newick_string(test_bad_ending))
# println(is_valid_newick_string(test_all_gut))
