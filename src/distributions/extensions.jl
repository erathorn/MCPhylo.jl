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


#################### Type Aliases ####################

const SymDistributionType = Union{
    Type{Biweight},
    Type{Cosine},
    Type{Epanechnikov},
    Type{Normal},
    Type{SymTriangularDist},
    Type{Triweight},
    Type{SymUniform},
}

const KernelDensityType = SymDistributionType
