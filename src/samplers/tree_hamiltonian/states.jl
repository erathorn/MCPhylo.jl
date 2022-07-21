#### Structs for positions in the hamiltonian space

abstract type HMC_State end

mutable struct Tree_HMC_State{T<:GeneralNode} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
    extended::Bool
    nni::Int
end

function Tree_HMC_State(tree::T, r::Vector{Float64}, g::Vector{Float64}, lf::Float64)::Tree_HMC_State{T} where T<:GeneralNode
    Tree_HMC_State{T}(tree, r, g, lf, false, 0)
end


mutable struct Array_HMC_State{T<:Array{<:Real}} <: HMC_State
    x::T
    r::Vector{Float64}
    g::Vector{Float64}
    lf::Float64
end

function transfer(s1::T)::T where T<:HMC_State
    T(deepcopy(s1.x), s1.r[:],s1.g[:],s1.lf, s1.extended, s1.nni)
end

function transfer(s1::T)::T where T<:Array_HMC_State
    T(s1.x[:], s1.r[:],s1.g[:],s1.lf)
end

hamiltonian(s::HMC_State)::Float64 = s.lf - 0.5 * turbo_dot(s.r, s.r)

function transfer!(s1::Tree_HMC_State{T}, s2::Tree_HMC_State{T}) where T
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


#### PNUTS Parameter Structs

struct NUTS_StepParams
    μ::Float64
    δ::Float64
    γ::Float64
    κ::Float64
    t0::Int
    τ::Float64
end

function update_step(p::N, μ::Float64)::N where N<:NUTS_StepParams
    NUTS_StepParams(μ, p.δ, p.γ, p.κ, p.t0, p.τ)
end

mutable struct NUTSstepadapter
    m::Int
    s_bar::Float64
    x_bar::Float64
    params::NUTS_StepParams
    metro_acc_prob::Float64
    avg_nni::Float64

    function NUTSstepadapter(m, s_bar, x_bar, μ, δ, γ, κ, t0, τ)
        new(m, s_bar, x_bar, NUTS_StepParams(μ, δ, γ, κ, t0, τ), 0.0, 0.0)
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
end

NUTSMeta() = NUTSMeta(0,0.0,0.0,0.0)

function update!(x::NUTSMeta, y::NUTSMeta)::Nothing
    x.nni += y.nni
    x.alpha += y.alpha
    x.nalpha += y.nalpha
    x.accnni += y.accnni
    nothing
end
