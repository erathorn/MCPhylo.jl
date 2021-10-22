#################### PDMatDistribution ####################

const PDMatDistribution = Union{InverseWishart, Wishart}

function unlist(d::PDMatDistribution, X::AbstractArray)
  n = dim(d)
  #y = [triu(trues(n,n))]
  #y = similar(X, Int(n * (n + 1) / 2))
  y = [X[i,j] for i in 1:n for j in 1:n if i>=j]
  #y = vec(X)
  # k = 0
  # for i in 1:n, j in i:n
  #   k += 1
  #   y[k] = X[i, j]
  # end
  y
end

function getix(row, col)::Int
  #assume left/lower triangular
  if col > row
      row = col
      col = row
  end
  (row-1)*row/2 + col
end


function relistlength(d::PDMatDistribution, X::AbstractArray{T}) where T
  n = dim(d)
  
  #Y = similar(X, n, n)
  #Y = zeros(n,n)
  #f(x,y) = Int(x * (y+1) >>> 2)
  #[i<=j ? X[f(i,j)] : 0 for i=1:n, j=1:n ]
  Y = [i>=j ? X[getix(i,j)] : X[getix(j,i)] for j=1:n, i=1:n]
  k = Int(n * (n + 1) / 2)
  (Y, k)
  #(collect(reshape(X, (n,n))), n*n)
end

# function link(d::PDMatDistribution, X::DenseMatrix)
#   n = dim(d)
#   Y = zeros(n, n)
#   U = cholesky(X).U
#   for i in 1:n
#     Y[i, i] = log(U[i, i])
#   end
#   for i in 1:n, j in (i + 1):n
#     Y[i, j] = U[i, j]
#   end
#   Y
# end

# function invlink(d::PDMatDistribution, X::DenseMatrix)
#   n = dim(d)
#   U = zeros(n, n)
#   di = diagind(n,n)
#   XD = diagm(exp.(X[di]) .- X[di])
#   U = X .+ XD
#   transpose(U) * U

#   # for i in 1:n
#   #   U[i, i] = exp(X[i, i])
#   # end
#   # for i in 1:n, j in (i + 1):n
#   #   U[i, j] = U[i,j] + X[i, j]
#   # end
#   # transpose(U) * U
# end

function logpdf(d::PDMatDistribution, X::DenseMatrix, transform::Bool)
  lp = logpdf(d, X)
  if transform && isfinite(lp)
    U = cholesky(X).U
    n = dim(d)
    for i in 1:n
      lp += (n - i + 2) * log(U[i, i])
    end
    lp += n * log(2)
  end
  lp
end
