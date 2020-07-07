
function pointwise_mat(arr1::Array{Float64,2},arr2::Array{Float64,2}, n_c::Int64)::Array{Float64,2}
    out::Array{Float64,2} = zeros(size(arr1))
    @simd for i in 1:n_c
        @inbounds out[1,i] = arr1[1,i] * arr2[1,i]
        @inbounds out[2,i] = arr1[2,i] * arr2[2,i]
    end # for
    return out
end # function pointwise_mat


function pointwise_mat!(arr1::Array{Float64,2},arr2::Array{Float64,2}, n_c::Int64)
    @simd for i in 1:n_c
        @inbounds arr1[1,i] = arr1[1,i] * arr2[1,i]
        @inbounds arr1[2,i] = arr1[2,i] * arr2[2,i]
    end # for
end # function pointwise_mat

function pointwise_vec(arr1::Array{Float64,1},arr2::Array{Float64,1})::Array{Float64}
    out::Array{Float64,1} = zeros(size(arr1))
    @simd for i = eachindex(out)
        @inbounds out[i] = arr1[i] * arr2[i]
    end # for
    return out
end # function pointwise_vec


function pointwise_vec(arr1::Array{Float64,2}, vec1::Vector{Float64}, vec2::Vector{Float64}, out::Vector{Float64})::Nothing
    @simd for i = eachindex(out)
        @inbounds out[i] = arr1[1, i]*vec1[i] + arr1[2, i]*vec2[i]
    end
    nothing
end # function


function my_dot(arr1::Array{Float64,2}, arr2::Array{Float64,1}, out_arr::Vector{Float64})::Nothing
    @simd for i = eachindex(out_arr)
        @inbounds out_arr[i] = arr1[1, i] * arr2[1] + arr1[2,i] * arr2[2]
    end #for
    nothing
end # function my_dot

function my_dot(arr1::Array{Float64,2}, arr2::Array{Float64,1}, n_c::Int64)::Array{Float64}
    out_arr = Array{Float64,1}(undef, n_c)
    @simd for i = eachindex(out_arr)
        @inbounds out_arr[i] = arr1[1, i] * arr2[1] + arr1[2,i] * arr2[2]
    end #for
    out_arr
end # function my_dot
