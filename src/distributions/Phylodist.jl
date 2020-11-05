
mutable struct PhyloDist <: DiscreteMatrixDistribution
    my_tree::T where T <: GeneralNode
    mypi::Float64
    rates::Vector{Float64}
    nbase::Int64
    nsites::Int64
    nnodes::Int64
    blv::Vector{Float64}

    function PhyloDist(my_tree::T, mypi::Float64, rates::Vector{Float64}, nbase::Int64, nsites::Int64, nnodes::Int64) where T <: GeneralNode
        new(my_tree, mypi, rates, nbase, nsites, nnodes, zeros(nnodes-1))
    end
end

function PhyloDist(my_tree::T, mypi::S, rates::A, nbase::Int64, nsites::Int64, nnodes::Int64) where {T<:TreeVariate, S<:ScalarVariate, A<:ArrayVariate}
    PhyloDist(my_tree.value, mypi.value, rates.value, nbase, nsites, nnodes)
end
minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nbase, d.nsites, d.nnodes)

function logpdf(d::PhyloDist, x::AbstractArray)
    mt = post_order(d.my_tree)

    blv = get_branchlength_vector(d.my_tree)
    return FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, blv)
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    blv = get_branchlength_vector(d.my_tree)
    mt = post_order(d.my_tree)

    f(y) = FelsensteinFunction(mt, d.mypi, d.rates, x, d.nsites, y)
    r = Zygote.pullback(f, blv)

    r[1], r[2](1.0)[1]
end
