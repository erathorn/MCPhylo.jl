module Converter

include("./Tree_Basics.jl")
using DataFrames
using ..Tree_Basics: Node, post_order

export to_df

function to_df(root::Node)#::DataFrame{Float64}

    post_order_iteration = post_order(root)

    name_list = [i.name for i in post_order_iteration]

    temp_ar = zeros(Float64, (length(post_order_iteration), length(post_order_iteration)))

    for i in post_order_iteration
        if i.nchild != 0
            ind = indexin(i.name, name_list)
            for j in i.child
                ind2 = indexin(j.name, name_list)
                temp_ar[ind[1], ind2[1]] = j.inc_length
            end
        end
    end
    df = convert(DataFrame, temp_ar)
    names!(df, [Symbol(i) for i in name_list])
    return df
end

end # module Converter
