
mutable struct PhyloDist <: DiscreteMatrixDistribution
    my_tree::T where T <: GeneralNode
    mypi::Array{Float64}
    rates::Array{Float64}
    nbase::Int64
    nsites::Int64
    nnodes::Int64
    blv::Vector{Float64}

    function PhyloDist(my_tree::T, mypi::Array{Float64}, rates::Vector{R}, nbase::Int64, nsites::Int64, nnodes::Int64) where {T <: GeneralNode, R<: Real}
        new(my_tree, mypi, rates, nbase, nsites, nnodes, zeros(nnodes-1))
    end

end

function PhyloDist(my_tree::T, mypi::S, rates::A, nbase::Int64, nsites::Int64, nnodes::Int64) where {T<:TreeVariate, S<:DenseArray{Float64}, A<:DenseArray{Float64}}
    PhyloDist(my_tree.value, Array(mypi), Array(rates), nbase, nsites, nnodes)
end

function PhyloDist(my_tree::T, mypi::S, rates::R, nbase::Int64, nsites::Int64, nnodes::Int64) where {T<:TreeVariate, S<:DenseArray{Float64}, R<:Real}
    PhyloDist(my_tree.value, Array(mypi), [rates], nbase, nsites, nnodes)
end


minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nbase, d.nsites, d.nnodes)

function logpdf(d::PhyloDist, x::AbstractArray)
    mt = post_order(d.my_tree)
    data = Array{Float64, 4}(undef, d.nbase, d.nsites, length(d.rates), d.nnodes)
    @inbounds for i in 1:length(d.rates)
        data[:, :, i, :] .= x
    end

    Felsenstein_grad(mt, d.mypi, d.rates, data, false)[1]
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.my_tree)
    data = Array{Float64, 4}(undef, d.nbase, d.nsites, length(d.rates), d.nnodes)
    @inbounds for i in 1:length(d.rates)
        data[:, :, i, :] .= x
    end

    resultate = Felsenstein_grad(mt, d.mypi, d.rates, data)
    #println(resultate)
    resultate
end
