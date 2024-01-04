using Distributions, Statistics, LinearAlgebra, LowRankEnKF

import Base: size, length

import Statistics: mean, var, std

#=
export initialize, initialize!

function initialize(N::Int, Distx0::MultivariateDistribution)
    NS = length(Distx0)

    # pre-allocate space
    X0 = 
    ENS = BasicEnsembleMatrix()
    ENS = EnsembleState(N, zeros(NS))

    ENS.S .= [rand(Distx0) for i = 1:N]
    return ENS
end
#
function initialize(N::Int, NS::Int)
    # pre-allocate space
    ENS = EnsembleState(N, zeros(NS))

    Dist = MvNormal(zeros(NS), I)
    ENS.S .= [rand(Dist) for i = 1:N]
    return ENS
end
#
function initialize!(X::BasicEnsembleMatrix{Nx, Ne, T})  where {Nx, Ne, T}
    Dist = MvNormal(zeros(Nx), I)
    for j in 1:Ne
        X(j) .= initialize(X(j))
    end
    return X
end