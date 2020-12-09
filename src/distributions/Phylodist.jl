
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

    blv = get_branchlength_vector(d.my_tree)
    res = Threads.Atomic{Float64}(0.0)

    Base.Threads.@threads for r in d.rates
        Threads.atomic_add!(
        res, FelsensteinFunction(mt, d.mypi, r, x, d.nsites, blv)
        )
    end
    res[]
end


function gradlogpdf(d::PhyloDist, x::AbstractArray)

    mt = post_order(d.my_tree)
    N = length(mt)-1
    mg = Threads.Atomic{Float64}.(zeros(N))
    res = Threads.Atomic{Float64}(0.0)


    # use ∇(F+G) = ∇F + ∇G to speed up the process
    Base.Threads.@threads for mrx in d.rates
        r1, gv = Felsenstein_grad(mt, d.mypi, mrx, x, d.nsites)
        Threads.atomic_add!(res, r1)
        @inbounds @simd for i in 1:N
            Threads.atomic_add!(mg[i], gv[i])
        end
    end

    res[], getproperty.(mg, :value)
end
