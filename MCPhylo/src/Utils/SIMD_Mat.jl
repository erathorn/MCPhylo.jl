
function pointwise_mat(arr1::Array{Float64,2},arr2::Array{Float64,2}, n_c::Int64)::Array{Float64}
    out::Array{Float64,2} = zeros(size(arr1))
    @inbounds for i in 1:n_c
        out[1,i] = arr1[1,i] * arr2[1,i]
        out[2,i] = arr1[2,i] * arr2[2,i]
    end # for
    return out
end # function pointwise_mat

function pointwise_vec(arr1::Array{Float64,1},arr2::Array{Float64,1}, n_c::Int64)::Array{Float64}
    out::Array{Float64,1} = zeros(size(arr1))
    @inbounds for i in 1:n_c
        out[i] = arr1[i] * arr2[i]
    end # for
    return out
end # function pointwise_vec



function my_dot(arr1::Array{Float64,2}, arr2::Array{Float64,1}, n_c::Int64)::Array{Float64}
    out::Array{Float64,1} = zeros(n_c)
    @inbounds for i in 1:n_c
        out[i] = arr1[1, i] * arr2[1] + arr1[2,i] * arr2[2]
    end #for
    return out
end # function my_dot
