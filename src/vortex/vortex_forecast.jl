export forecast, VortexForecast




#### FORECAST OPERATORS ####


mutable struct VortexForecast{Nx,withfreestream,BT} <: AbstractForecastOperator{Nx}
    vm :: VortexModel
    pfb :: PotentialFlowBody
end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nx = 3*length(vm.vortices)
    body = pfb.points
    VortexForecast{Nx,withfreestream,typeof(body)}(vm,pfb)
end


function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{Nx}) where {Nx}
    @unpack vm, pfb = fdata
    @unpack points = pfb
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices[end-1:end])
    pushvortices!(vm,vLEnew,vTEnew)
    
    xnew = x[1:end]
    append!(xnew,zeros(6))
    for (i, vortex) in enumerate(vm.vortices)
        xnew[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
    
    return xnew
end


# vortices released at one-third of the way from the edge to the last released vortex from that edge
function createsheddedvortices(plate::Plate,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE, vTE
end