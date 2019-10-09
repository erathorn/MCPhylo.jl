
mutable struct PhyloDist <: DiscreteMatrixDistribution
    my_tree
    mypi::Real
    #rates::Vector
    nnodes::Int64
    nbase::Int64
    nsites::Int64
end
minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

Base.size(d::PhyloDist) = (d.nnodes, d.nbase, d.nsites)

function logpdf(d::PhyloDist, x::AbstractArray)
    
    mt = post_order(d.my_tree.value)
    rates = ones(3132)
    return FelsensteinFunction(mt, d.mypi, rates, x, d.nsites)
end

function gradlogpdf(d::PhyloDist, x::AbstractArray)
    mt = pre_order(d.my_tree.value)
    return GradiantLog(mt, d.mypi, ones(3132), x, d.nsites)
end
