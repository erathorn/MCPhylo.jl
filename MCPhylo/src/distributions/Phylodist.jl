
mutable struct PhyloDist <: ContinuousUnivariateDistribution
    my_tree
    mypi::Real
    rates::Vector
end
minimum(d::PhyloDist) = -Inf
maximum(d::PhyloDist) = Inf

function logpdf(d::PhyloDist, x::Real)

    mt = post_order(d.my_tree.value)
    return FelsensteinFunction(mt, d.mypi, d.rates)
end

function gradlogpdf(d::PhyloDist, x::Real)
    mt = pre_order(d.my_tree.value)
    return GradiantLog(mt, d.mypi, d.rates)
end
