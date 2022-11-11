#################### Flat Distribution ####################

struct Flat <: ContinuousUnivariateDistribution end

minimum(d::Flat) = -Inf
maximum(d::Flat) = Inf
insupport(d::Flat, x::Real) = true

logpdf(d::Flat, x::Real) = 0.0




#################### Null Distribution ####################

struct NullUnivariateDistribution <: UnivariateDistribution{ValueSupport} end


#################### SymUniform Distribution ####################

struct SymUniform <: ContinuousUnivariateDistribution
    SymUniform() = Uniform(-1.0, 1.0)
    SymUniform(μ::Real, σ::Real) = Uniform(μ - σ, μ + σ)
end

struct SymBinary{R<:Real} <: DiscreteUnivariateDistribution
    w::R

    SymBinary() = new{Float64}(0.5)
    SymBinary(w::F) where {F<:Real} = new{F}(w)
    SymBinary(x, y) = SymBinary()
end

function rand(r::R, d::SymBinary{F}) where {R<:AbstractRNG,F<:Real}
    rand(r) > d.w ? 1 : -1
end

Base.maximum(d::SymBinary{F}) where {F} = 1
Base.minimum(d::SymBinary{F}) where {F} = -1

#################### Type Aliases ####################

const SymDistributionType = Union{
    Type{Biweight},
    Type{Cosine},
    Type{Epanechnikov},
    Type{Normal},
    Type{SymTriangularDist},
    Type{Triweight},
    Type{SymUniform},
    Type{SymBinary},
}

const KernelDensityType = SymDistributionType
