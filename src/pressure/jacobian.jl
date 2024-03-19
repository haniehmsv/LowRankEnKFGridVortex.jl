export jacob!

function jacob!(H,x::AbstractVector,t,Δt,obs::AbstractObservationOperator{Ny}) where {Ny}
    dhdx(x,t,Δt,obs) = ForwardDiff.jacobian(v->observations(v,t,Δt,obs),x)
    H .= dhdx(x,t,Δt,obs)
end