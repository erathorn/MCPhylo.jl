#### Structs for positions in the hamiltonian space

abstract type HMC_State end

mutable struct Tree_HMC_State{T<:GeneralNode} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
    extended::Bool
    nni::Int
    att_nni::Int
end



function Tree_HMC_State(
    tree::T,
    r::Vector{Float64},
    g::Vector{Float64},
    lf::Float64,
)::Tree_HMC_State{T} where {T<:GeneralNode}
    Tree_HMC_State{T}(tree, r, g, lf, false, 0, 0)
end


struct Extendend_Tree_HMC_State{H<:Tree_HMC_State}
    curr_state::H
    ext_state::H
    extended::Vector{Bool}
    direction::Int
end

function Extendend_Tree_HMC_State(tree, r, g, lf, direction)
    t = Tree_HMC_State(deepcopy(tree), deepcopy(r), deepcopy(g), lf)
    t2 = Tree_HMC_State(deepcopy(tree), deepcopy(r), deepcopy(g), lf)
    Extendend_Tree_HMC_State(t, t2, [false], direction)
end


mutable struct Array_HMC_State{T<:Array{<:Real}} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
end

function transfer(s1::T)::T where {T<:HMC_State}
    T(deepcopy(s1.x), s1.r[:], s1.g[:], s1.lf, s1.extended, s1.nni, s1.att_nni)
end

function transfer(s1::T)::T where {T<:Array_HMC_State}
    T(s1.x[:], s1.r[:], s1.g[:], s1.lf)
end

hamiltonian(s::HMC_State)::Float64 = s.lf - 0.5 * turbo_dot(s.r, s.r)

function transfer!(s1::Tree_HMC_State{T}, s2::Tree_HMC_State{T}) where {T}
    s1.x = deepcopy(s2.x)
    s1.r[:] .= s2.r[:]
    s1.g[:] .= s2.g[:]
    s1.lf = s2.lf
    s1.extended = s2.extended
    s1.nni = s2.nni
    nothing
end

function transfer!(s1::Array_HMC_State, s2::Array_HMC_State)
    s1.x[:] .= s2.x[:]
    s1.r[:] .= s2.r[:]
    s1.g[:] .= s2.g[:]
    s1.lf = s2.lf
    nothing
end



function unextend!(S::Extendend_Tree_HMC_State)
    transfer!(S.curr_state, S.ext_state)
    S.extended[1] = false
    nothing
end

function extend!(S::Extendend_Tree_HMC_State{P}, PS::P) where {P}
    transfer!(S.ext_state, PS)
    S.extended[1] = true
    nothing
end

function transfer(S::Extendend_Tree_HMC_State)
    deepcopy(S)
end

function transfer!(S1::Extendend_Tree_HMC_State, S2::Extendend_Tree_HMC_State)
    @assert S1.direction == S2.direction
    transfer!(S1.curr_state, S2.curr_state)
    transfer!(S1.ext_state, S2.ext_state[1])
    S1.extended[1] = S2.extended[1]
    nothing
end


hamiltonian(S::Extendend_Tree_HMC_State) = hamiltonian(S.curr_state)



#### PNUTS Parameter Structs

struct NUTS_StepParams
    μ::Float64
    δ::Float64
    γ::Float64
    κ::Float64
    t0::Int
    δ_NNI::Float64
end

function update_step(p::N, μ::Float64)::N where {N<:NUTS_StepParams}
    NUTS_StepParams(μ, p.δ, p.γ, p.κ, p.t0, p.δ_NNI)
end

mutable struct NUTSstepadapter
    m::Int
    s_bar::Float64
    x_bar::Float64
    params::NUTS_StepParams
    metro_acc_prob::Float64
    NNI_stat::Float64

    function NUTSstepadapter(m, s_bar, x_bar, μ, δ, δ_NNI, γ, κ, t0)
        new(m, s_bar, x_bar, NUTS_StepParams(μ, δ, γ, κ, t0, δ_NNI), 0.0, 0.0)
    end

    function NUTSstepadapter(m, s_bar, x_bar, params::NUTS_StepParams)
        new(m, s_bar, x_bar, params, 0.0, 0.0)
    end
end


mutable struct NUTSMeta
    nni::Int
    alpha::Float64
    nalpha::Float64
    accnni::Int
    att_nni::Int
    l_NNI::Int
end

NUTSMeta() = NUTSMeta(0, 0.0, 0.0, 0, 0, 0)
