export jacob!

function jacob!(H,x::AbstractVector,t,Δt,obs::AbstractObservationOperator{Ny},i::Int) where {Ny}
    dhdx(x,t,Δt,obs,i) = ForwardDiff.jacobian(v->observations(v,t,Δt,obs,i),x)
    H .= dhdx(x,t,Δt,obs,i)
end

jacob!(H,x::AbstractVector,t,Δt,obs::AbstractObservationOperator{Ny}) where {Ny} = nothing