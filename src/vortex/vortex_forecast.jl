export forecast, VortexForecast

import LowRankEnKF: forecast


#### FORECAST OPERATORS ####


mutable struct VortexForecast{Nx,withfreestream,BT,Ne} <: AbstractForecastOperator{Nx}

    "vortex model from GridPotentialFlow.jl"
    vm :: VortexModel

    "potential flow body from GridPotentialFlow.jl"
    pfb :: PotentialFlowBody

    "Number of vortices"
    Nv :: Int64

end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vm::VortexModel{Nb,Ne},pfb::PotentialFlowBody) where {Nb,Ne}
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nv = length(vm.vortices)
    Nx = 3*Nv
    body = pfb.points
    VortexForecast{Nx,withfreestream,typeof(body),Ne}(vm,pfb,Nv)
end


function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{Nx,Ne}) where {Nx,Ne}
    @unpack vm, pfb, Nv = fdata
    @unpack points = pfb
    time_advancement!(vm,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices[end-1:end])
    pushvortices!(vm,vLEnew,vTEnew)

    xnew = x[1:end]
    # New vortices released from the two edges augment the state vector by 3*Ne
    append!(xnew,zeros(6))
    vortices_to_states!(xnew,vm)
    Nv = length(vm.vortices)
    return xnew
end

function time_advancement!(vm::VortexModel,Δt)
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

# vortices released at one-third of the way from the edge to the last released vortex from that edge
function createsheddedvortices(plate::Plate,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE, vTE
end
