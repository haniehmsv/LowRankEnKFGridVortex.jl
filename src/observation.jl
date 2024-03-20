# subtypes of AbstractObservationOperator should provide
#  - an extension of the observation operator
#  - an extension of the jacob! operator
# They should keep sensor locations and vortex configuration internal
# Nx is number of states, Ny is number of observations

export observations!, AbstractObservationOperator, jacob!,
       state_filter!, normal_loglikelihood, log_uniform, measurement_length, state_length

abstract type AbstractObservationOperator{Ny,withsensors} end


measurement_length(::AbstractObservationOperator{Ny}) where {Ny} = Ny
state_length(X::EnsembleMatrix) = size(X,1)
state_length(config::VortexForecast) = config.Nv*3


"""
    normal_loglikelihood(x,t,ystar,Σϵ,obs) -> Float64

For a given state `x` at time `t`, return the log of the likelihood function,
given observations `ystar`, noise covariance `Σϵ`, and observation structure `obs`.
"""
function normal_loglikelihood(x,t,ystar,Σϵ,obs::AbstractObservationOperator)
    y = observations(x,t,obs)
    #loss = norm(ystar-y .- mean(ystar-y),Σϵ)
    loss = norm(ystar-y,Σϵ)
    return -loss^2/2
end

"""
    log_uniform(x,bounds)

Given state vector `x` and bounds vector `bounds`, return
the log of the scaled uniform probability of `x` relative to the
bounds. It returns -Inf if `x` is outside the bounds and 0 if `x` is inside
or on the bounds.
"""
function log_uniform(x::Vector,bounds::Vector{<:Tuple})
    logp = 0.0
    for i in eachindex(x)
        a, b = bounds[i]
        logp += log((a<=x[i]<=b)*1.0)
    end
    return logp
end

# """
#     observations(x::AbstractVector,t::Float64,obs::AbstractObservationOperator) -> X

# Compute the observation function `h` for state `x` at time `t`.
# The function `h` should take as inputs a vector of measurement points (`sens`), a vector of vortices,
# and the configuration data `config`.
# """
# function observations(x::AbstractVector,t,Δt,obs::AbstractObservationOperator)
#     return observations(x,t,Δt,obs)
# end

"""
    observations!(Y::EnsembleMatrix,X::EnsembleMatrix,t::Float64,obs::AbstractObservationOperator)

Compute the observation function `h` for each of the states in `X` and place them in `Y`.
The function `h` should take as inputs a Ny-dimensional vector of measurement points (`sens`), a vector of vortices,
and the configuration data `config`.
"""
function observations!(Y::EnsembleMatrix{Ne},X::EnsembleMatrix{Ne},t,Δt,obs::AbstractObservationOperator{Ny}) where {Ny,Ne}
  for j in 1:Ne
      Y(j) .= observations(X(j),t,Δt,obs)
  end
  return Y
end

"""
    jacob!(J,x::AbstractVector,t::Float64,obs::AbstractObservationOperator)

Compute the Jacobian of the observation function at state `x` and return it in `J`.
"""
function jacob!(J,x::AbstractVector,t,obs::AbstractObservationOperator) end


state_filter!(x,obs::AbstractObservationOperator) = x
