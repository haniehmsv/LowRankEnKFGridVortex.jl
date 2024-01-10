import Base: size, length

import Statistics: mean, var, std

export initialize!

function initialize!(X::BasicEnsembleMatrix{Ne, T}; Σi=I)  where {Ne, T}
    Dist = [MvNormal(X.X[:,j], Σi) for j in 1:Ne]
    for j in 1:Ne
        X(j) .= rand(Dist[j])
    end
    return X
end