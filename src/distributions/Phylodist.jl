
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

    Nbases, Nsites, Nrates, Nnodes = size(data)
    Q::Array{Float64,2} = ones(Nbases,Nbases)
    Q[diagind(Nbases,Nbases)] .= -1
    Q .*= d.mypi
    D, _ , U, _ = LAPACK.geev!('V', 'N',Q)
    Uinv = inv(U)

    mu::Float64 =  1.0 / (2.0 * prod(d.mypi))
    mutationArray::Array{Float64,4} = Array{Float64,4}(undef, Nbases, Nbases, Nrates, Nnodes-1)


    Down::Array{Float64,4} = similar(data)
    fels_ll(mt, data, D, U, Uinv, d.rates, mu, Nrates, Nsites, Down, d.mypi, mutationArray)
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.my_tree)
    data = Array{Float64, 4}(undef, d.nbase, d.nsites, length(d.rates), d.nnodes)
    @inbounds for i in 1:length(d.rates)
        data[:, :, i, :] .= x
    end
    # use ∇(F+G) = ∇F + ∇G to speed up the process
    #Base.Threads.@threads for mrx in d.rates
    Felsenstein_grad(mt, d.mypi, d.rates, data)
    #    Threads.atomic_add!(res, r1)
    #    @inbounds @simd for i in 1:N
    #        Threads.atomic_add!(mg[i], gv[i])
    #    end
    #end

    #res[], getproperty.(mg, :value)
end
