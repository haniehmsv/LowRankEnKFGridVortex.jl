

export forecast, VortexForecast, advect_vortices!, createsheddedvortices, construct_intermediate_model!, retrieve_vm_from_intermediatevm!

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
function VortexForecast(vvm::Vector{VortexModel{Nb,Ne,TS,TU,TE,TF,TX,ILS}}) where {Nb,Ne,TS,TU,TE,TF,TX,ILS}
    withfreestream = vvm[1].U∞ == 0.0 ? false : true
    Nv = length(vvm[1].vortices)
    Nx = 3*Nv
    VortexForecast{withfreestream,Nb,Ne}(vvm,Nv)
end



function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{Nb,Ne},i::Int64) where {Nb,Ne}
    @unpack vvm = fdata
    vm = vvm[i] #i-th ensemble member
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb
    
    states_to_vortices!(vm,x,Δt)
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

"""construct_intermediate_model!(intermediate_vm::VortexModel{Nb,0},vm::VortexModel{Nb,Ne}) --> VortexModel{Nb,0}
an intermediate vortex model with all fields the same as vm except that intermediate_vm has no regularized edge. This allows
the solution of the system for the existing vortices which solves solve!(sol::ConstrainedIBPoissonSolution, 
vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}}) in the vortexmodel.jl file. 
"""
function construct_intermediate_model!(intermediate_vm::VortexModel{Nb,0},vm::VortexModel{Nb,Ne}) where {Nb,Ne}
    intermediate_vm.bodies = deepcopy(vm.bodies)
    intermediate_vm.vortices = deepcopy(vm.vortices)
    intermediate_vm.U∞ = deepcopy(vm.U∞)
    for i=1:Nb
        intermediate_vm.bodies[i].edges = Int64[]
    end
end

function retrieve_vm_from_intermediatevm!(vm::VortexModel{Nb,Ne},intermediate_vm::VortexModel{Nb,0}) where {Nb,Ne}
    getΓ.(vm.bodies) .= deepcopy(getΓ.(intermediate_vm.bodies))
    vm.vortices = deepcopy(intermediate_vm.vortices)
    vm.U∞ = deepcopy(intermediate_vm.U∞)
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain and a body with Ne=1 regularized edge.
Used in the foreward model."""
function advect_vortices!(vm::VortexModel{Nb,Ne},Δt) where {Nb,Ne}
    X = getvortexpositions(vm)
    subtractcirculation!(vm.bodies, vm.vortices.Γ[end-Ne])
    Ẋ = vortexvelocities!(vm)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain and a body with 0 regularized edge.
Used in the observation model."""
function advect_vortices!(vm::VortexModel{Nb,0},Δt)
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