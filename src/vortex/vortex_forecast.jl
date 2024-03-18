
import GridPotentialFlow: AbstractPotentialFlowSystem, UnsteadyRegularizedIBPoisson, ConstrainedIBPoisson
export forecast, VortexForecast, advect_vortices!, createsheddedvortices

#### FORECAST OPERATORS ####


mutable struct VortexForecast{withfreestream,Nb,Ne,TS<:Union{AbstractPotentialFlowSystem,Laplacian}} <: AbstractForecastOperator

    "vortex model from GridPotentialFlow.jl"
    vvm :: Vector{<:VortexModel{Nb,Ne,TS}}

    "Number of vortices"
    Nv :: Int64

end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vvm::Vector{<:VortexModel{Nb,Ne,TS}}) where {Nb,Ne,TS}
    withfreestream = vvm[1].U∞ == 0.0 ? false : true
    Nv = length(vvm[1].vortices)
    Nx = 3*Nv
    VortexForecast{withfreestream,Nb,Ne,TS}(vvm,Nv)
end


"""System with regularized edges. Enforce circulation constraints."""
function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{true,Nb,Ne,<:UnsteadyRegularizedIBPoisson{Nb,Ne}},i::Int64) where {Nb,Ne}
    @unpack vvm = fdata
    vm = vvm[i] #i-th ensemble member
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb
    
    states_to_vortices!(vm,x,Δt)
    vm.bodies[1].Γ = -sum(vm.vortices.Γ[1:end-1])
    advect_vortices!(vm,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices,DT=Real)
    pushvortices!(vm,vLEnew,vTEnew)
    xnew = similar(x[1:end])
    # New vortices released from the two edges augment the state vector by 3*Ne
    append!(xnew,zeros(6))
    vortices_to_states!(xnew,vm,x)
    fdata.Nv = length(vm.vortices)
    return xnew
end

"""System without regularized edges. Enforce circulation constraints."""
function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{true,Nb,Ne,<:ConstrainedIBPoisson{Nb}},i::Int64) where {Nb,Ne}
    @unpack vvm = fdata
    vm = vvm[i] #i-th ensemble member
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    states_to_vortices!(vm,x)
    vm.bodies[1].Γ = -sum(vm.vortices.Γ)
    advect_vortices!(vm,Δt)
    vortices_to_states!(x,vm)
    return x
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
function createsheddedvortices(plate::Polygon,oldvortices;DT=Float64)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0,DT=DT)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0,DT=DT)

    return vLE, vTE
end