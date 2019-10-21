function setinits!(d::TreeStochastic, m::Model, x::Node)
    d.value = x
    d.distr = d.eval(m)
    insupport(d.distr, x) || throw(ArgumentError("The supplied tree does not match the topological tree constraints."))
    setmonitor!(d, d.monitor)
end # function
