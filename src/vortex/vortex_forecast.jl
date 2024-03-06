

export forecast, VortexForecast, advect_vortices!, createsheddedvortices

#### FORECAST OPERATORS ####


mutable struct VortexForecast{withfreestream,Nb,Ne} <: AbstractForecastOperator

    "vortex model from GridPotentialFlow.jl"
    vvm :: Vector{VortexModel}

    "Number of vortices"
    Nv :: Int64

end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vvm::Vector{<:VortexModel{Nb,Ne}}) where {Nb,Ne}
    withfreestream = vvm[1].U∞ == 0.0 ? false : true
    Nv = length(vvm[1].vortices)
    Nx = 3*Nv
    VortexForecast{withfreestream,Nb,Ne}(vvm,Nv)
end


function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{withfreestream,Nb,Ne},i::Int64) where {withfreestream,Nb,Ne}
    @unpack vvm = fdata
    vm = vvm[i] #i-th ensemble member
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb
    
    states_to_vortices!(vm,x,Δt)
    vm.bodies[1].Γ = -sum(vm.vortices.Γ[1:end-1])
    advect_vortices!(vm,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices[end-1:end])
    pushvortices!(vm,vLEnew,vTEnew)
    xnew = similar(x[1:end])
    # New vortices released from the two edges augment the state vector by 3*Ne
    append!(xnew,zeros(6))
    vortices_to_states!(xnew,vm,x,Δt)
    fdata.Nv = length(vm.vortices)
    return xnew
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain and a body with Ne=1 regularized edge. Also solves for the new vortex shedded at the TE.
Used in the foreward model."""
function advect_vortices!(vm::VortexModel{Nb,Ne},Δt) where {Nb,Ne}
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain.
Used in the observation model."""
function advect_vortices!(vm::VortexModel{Nb,Ne},sol::ConstrainedIBPoissonSolution,Δt) where {Nb,Ne}
    X = getvortexpositions(vm)
    Ẋ = deepcopy(X)
    vortexvelocities!(Ẋ, vm, sol.ψ)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

# vortices released at one-third of the way from the edge to the last released vortex from that edge
function createsheddedvortices(plate::Polygon,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*FD.value(oldvortices[end-1].x),2/3*plate.y[1]+1/3*FD.value(oldvortices[end-1].y),0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*FD.value(oldvortices[end].x),2/3*plate.y[end]+1/3*FD.value(oldvortices[end].y),0.0)

    return vLE, vTE
end