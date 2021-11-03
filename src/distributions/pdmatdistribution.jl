#################### PDMatDistribution ####################

function unlist(d::PDMatDistribution, X::AbstractArray)
    n = dim(d)
    y = [X[i, j] for i = 1:n for j in 1:n if i >= j]
    y
end

function getix(row, col)::Int
    #assume left/lower triangular
    if col > row
        row = col
        col = row
    end
    (row - 1) * row / 2 + col
end


function relistlength(d::PDMatDistribution, X::AbstractArray{T}) where {T}
    n = dim(d)
    Y = [i >= j ? X[getix(i, j)] : X[getix(j, i)] for j = 1:n, i = 1:n]
    k = Int(n * (n + 1) / 2)
    (Y, k)
end
